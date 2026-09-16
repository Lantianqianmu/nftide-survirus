#!/usr/bin/env python
"""Auditable filtering of SurVirus results.remapped.txt (Python 2.7/3 compatible)."""
from __future__ import print_function
import argparse
import collections
import ctypes
import json
import math
import os
import re
import subprocess
import tempfile


def fasta(path):
    records = collections.OrderedDict()
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError('Duplicate FASTA record: ' + name)
                records[name] = ''
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence without FASTA header')
                records[name] += line.strip().upper()
    by_id = {}
    for name, seq in records.items():
        call_id = int(name.rsplit('_', 1)[0])
        if call_id in by_id:
            raise ValueError('Duplicate ID: ' + name)
        by_id[call_id] = (name, seq)
    return by_id


def masked_fractions(path, sequences):
    """Union of BED intervals, using standard zero-based, half-open coordinates."""
    intervals = collections.defaultdict(list)
    lengths = dict((name, len(seq)) for name, seq in sequences.values())
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            name, start, end = line.split()[:3]
            start, end = int(start), int(end)
            if name not in lengths or not 0 <= start <= end <= lengths[name]:
                raise ValueError('Invalid masking BED interval: ' + line.strip())
            intervals[name].append((start, end))
    result = {}
    for call_id, (name, seq) in sequences.items():
        covered, right = 0, 0
        for start, end in sorted(intervals[name]):
            covered += max(0, end - max(start, right))
            right = max(right, end)
        result[call_id] = covered / float(len(seq)) if seq else 0.0
    return result


def breakpoint(value):
    chrom, strand, start, end = value.rsplit(':', 3)
    if strand not in ('+', '-'):
        raise ValueError('Invalid breakpoint: ' + value)
    return chrom, strand, int(start) if strand == '-' else int(end)


def sdust_fractions(sequences, executable):
    """Mask ACGT runs independently: sdust 0.1 mis-offsets ends after N.

    Ambiguous bases break runs and remain unmasked; the denominator is the
    complete sequence length. Never clamp corrupt legacy BED coordinates.
    """
    runs = {}
    covered = collections.defaultdict(int)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa') as handle:
        for call_id, (name, seq) in sorted(sequences.items()):
            for match in re.finditer('[ACGT]+', seq):
                key = str(len(runs))
                runs[key] = (call_id, match.group())
                handle.write('>%s\n%s\n' % (key, match.group()))
        handle.flush()
        if runs:
            output = subprocess.check_output([executable, handle.name]).decode('ascii')
            intervals = collections.defaultdict(list)
            for line in output.splitlines():
                key, start, end = line.split()[:3]
                start, end = int(start), int(end)
                if key not in runs or not 0 <= start <= end <= len(runs[key][1]):
                    raise ValueError('Invalid recomputed sdust interval: ' + line)
                intervals[key].append((start, end))
            for key, spans in intervals.items():
                right = 0
                for start, end in sorted(spans):
                    covered[runs[key][0]] += max(0, end - max(start, right))
                    right = max(right, end)
    return dict((i, covered[i] / float(len(seq)) if seq else 0.0)
                for i, (name, seq) in sequences.items())


def pair_distance(a, b, minimum, maximum):
    ah, av, bh, bv = a['hbp'], a['vbp'], b['hbp'], b['vbp']
    if ah[0] != bh[0] or av[0] != bv[0] or ah[1] == bh[1] or av[1] == bv[1]:
        return None
    distance = ah[2] - bh[2] if ah[1] == '-' else bh[2] - ah[2]
    return abs(distance) if minimum <= distance <= maximum else None


def reverse_complement(seq):
    return ''.join(dict(zip('ACGTN', 'TGCAN')).get(c, 'N') for c in seq[::-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--ssw-library')
    parser.add_argument('--sdust', help='Recompute masks safely from ACGT runs using this executable')
    parser.add_argument('--min-host-pbs', type=float, default=0.8)
    parser.add_argument('--min-pairs', type=int, default=2)
    parser.add_argument('--split-min-pairs', type=int, default=1)
    parser.add_argument('--min-split-reads', type=int, default=1)
    parser.add_argument('--min-seq-length', type=int, default=30)
    parser.add_argument('--max-masked-fraction', type=float, default=0.8)
    parser.add_argument('--mask-mode', choices=['reject', 'flag'], default='reject')
    parser.add_argument('--pairing-rescue', choices=['true', 'false'], default='true')
    parser.add_argument('--pair-min-distance', type=int, default=-50)
    parser.add_argument('--pair-max-distance', type=int, default=1000)
    parser.add_argument('--deduplicate', choices=['true', 'false'], default='true')
    parser.add_argument('--dedup-fraction', type=float, default=0.8)
    args = parser.parse_args()
    for value in (args.min_host_pbs, args.max_masked_fraction, args.dedup_fraction):
        if not 0 <= value <= 1:
            parser.error('Fraction thresholds must be between 0 and 1')
    if min(args.min_pairs, args.split_min_pairs, args.min_split_reads, args.min_seq_length) < 1:
        parser.error('Support and sequence thresholds must be positive')
    if args.pair_min_distance > args.pair_max_distance:
        parser.error('Pair distance minimum exceeds maximum')
    if args.deduplicate == 'true' and not args.ssw_library:
        parser.error('--ssw-library is required when deduplication is enabled')
    h = fasta(os.path.join(args.workdir, 'host_bp_seqs.fa'))
    v = fasta(os.path.join(args.workdir, 'virus_bp_seqs.fa'))
    if args.sdust:
        hm = sdust_fractions(h, args.sdust)
        vm = sdust_fractions(v, args.sdust)
    else:
        hm = masked_fractions(os.path.join(args.workdir, 'host_bp_seqs.masked.bed'), h)
        vm = masked_fractions(os.path.join(args.workdir, 'virus_bp_seqs.masked.bed'), v)
    calls, seen = [], set()
    with open(os.path.join(args.workdir, 'results.remapped.txt')) as handle:
        for line in handle:
            if not line.strip():
                continue
            x = line.split()
            if len(x) != 13:
                raise ValueError('Expected 13 candidate fields: ' + line.strip())
            i, pairs, splits = int(x[0]), int(x[4]), int(x[5])
            if min(i, pairs, splits) < 0:
                raise ValueError('Negative ID/support count: ' + line.strip())
            pbs = float(x[7])
            if i in seen or i not in h or i not in v:
                raise ValueError('Duplicate ID or missing sequence: ' + x[0])
            seen.add(i)
            reasons = []
            if math.isnan(pbs) or math.isinf(pbs):
                reasons.append('INVALID_HOST_PBS')
            elif pbs < args.min_host_pbs:
                reasons.append('LOW_HOST_PBS')
            if not (pairs >= args.min_pairs or (pairs >= args.split_min_pairs and splits >= args.min_split_reads)):
                reasons.append('INSUFFICIENT_SUPPORT')
            if len(h[i][1]) < args.min_seq_length:
                reasons.append('SHORT_HOST_SEQUENCE')
            if len(v[i][1]) < args.min_seq_length:
                reasons.append('SHORT_VIRUS_SEQUENCE')
            repeat = []
            if hm[i] > args.max_masked_fraction:
                repeat.append('HOST_LOW_COMPLEXITY')
            if vm[i] > args.max_masked_fraction:
                repeat.append('VIRUS_LOW_COMPLEXITY')
            if args.mask_mode == 'reject':
                reasons.extend(repeat)
            calls.append(dict(id=i, raw=x, hbp=breakpoint(x[1]), vbp=breakpoint(x[2]),
                              reasons=reasons, repeat=repeat, status='rejected' if reasons else 'accepted',
                              partner='', duplicate='', rescued_by=''))

    # Pair strong calls first. Only original passing calls can rescue one
    # rejected partner; no chains of rescue. Ties follow candidate file order.
    initial = [c for c in calls if c['status'] == 'accepted']
    rejected = [c for c in calls if c['status'] == 'rejected']
    if args.pairing_rescue == 'true':
        for pool in (initial, rejected):
            for a in initial:
                if a['partner'] != '':
                    continue
                eligible = [(pair_distance(a, b, args.pair_min_distance, args.pair_max_distance), n, b)
                            for n, b in enumerate(pool) if b is not a and b['partner'] == '']
                eligible = [e for e in eligible if e[0] is not None]
                if eligible:
                    b = min(eligible, key=lambda e: (e[0], e[1]))[2]
                    a['partner'], b['partner'] = b['id'], a['id']
                    if b['status'] == 'rejected':
                        b['status'] = 'rescued'
                        b['rescued_by'] = a['id']

    if args.deduplicate == 'true':
        lib = ctypes.CDLL(os.path.abspath(args.ssw_library))
        lib.aligned_fraction.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        lib.aligned_fraction.restype = ctypes.c_double
    def similar(a, b, sequences, side):
        first, second = sequences[a['id']][1], sequences[b['id']][1]
        ref, query = (first, second) if len(first) >= len(second) else (second, first)
        if a[side][1] != b[side][1]:
            query = reverse_complement(query)
        return lib.aligned_fraction(ref.encode('ascii'), query.encode('ascii')) >= args.dedup_fraction
    kept = []
    for call in calls:
        if call['status'] == 'rejected':
            continue
        if args.deduplicate == 'true':
            for previous in kept:
                if similar(previous, call, h, 'hbp') and similar(previous, call, v, 'vbp'):
                    call['status'], call['duplicate'] = 'rejected', previous['id']
                    break
        if call['duplicate'] == '':
            kept.append(call)

    # A rescued call must retain its actual qualifying sponsor after dedup.
    # Keep the original rescue link for auditing, but clear stale pair links.
    kept_ids = set(c['id'] for c in kept)
    for call in kept:
        if call['status'] == 'rescued' and call['rescued_by'] not in kept_ids:
            call['status'] = 'rejected'
            call['reasons'].append('RESCUE_PARTNER_REMOVED')
    kept = [c for c in kept if c['status'] != 'rejected']
    kept_ids = set(c['id'] for c in kept)
    for call in calls:
        if call['id'] not in kept_ids or call['partner'] not in kept_ids:
            call['partner'] = ''

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    header = ['ID','HOST_BREAKPOINT','VIRUS_BREAKPOINT','SUPPORTING_PAIRS','SPLIT_READS','HOST_PBS',
              'COVERAGE','HOST_SEQ_LENGTH','VIRUS_SEQ_LENGTH','HOST_MASKED_FRACTION','VIRUS_MASKED_FRACTION',
              'REPEAT_FLAGS','STATUS','REJECTION_REASON','INITIAL_FILTER_FAILURES','PAIRED_WITH','DUPLICATE_OF',
              'RESCUED_BY','HOST_HAS_MASKED_BASES','VIRUS_HAS_MASKED_BASES']
    handles = dict((key, open(os.path.join(args.outdir, key + '.tsv'), 'w'))
                   for key in ('all_candidates', 'accepted', 'rejected'))
    for handle in handles.values():
        handle.write('\t'.join(header) + '\n')
    for call in calls:
        i, x = call['id'], call['raw']
        reason = ('DUPLICATE_SEQUENCE' if call['duplicate'] != '' else ';'.join(call['reasons'])) if call['status'] == 'rejected' else ''
        row = [i, x[1], x[2], x[4], x[5], x[7], (float(x[11])+float(x[12]))/2,
               len(h[i][1]), len(v[i][1]), hm[i], vm[i], ';'.join(call['repeat']), call['status'],
               reason, ';'.join(r for r in call['reasons'] if r != 'RESCUE_PARTNER_REMOVED'),
               call['partner'], call['duplicate'], call['rescued_by'], int(hm[i] > 0), int(vm[i] > 0)]
        text = '\t'.join(str(z) for z in row) + '\n'
        handles['all_candidates'].write(text)
        handles['rejected' if call['status'] == 'rejected' else 'accepted'].write(text)
    for handle in handles.values():
        handle.close()
    summary = dict(total=len(calls), initially_accepted=len(initial), retained=len(kept),
                   statuses=dict(collections.Counter(c['status'] for c in calls)), parameters=vars(args))
    with open(os.path.join(args.outdir, 'summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
