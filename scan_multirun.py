#!/usr/bin/env python3
"""Scan a BIDS dataset for sessions that hold several runs of the same scan.

Called by scan_multirun.sh, which passes the settings from pipeline_config.cfg.
Writes run_selection.tsv: one row per subject / session / scan type that has
more than one run, with a suggestion you can correct.

WHY A SUGGESTION AND NOT JUST A LIST
------------------------------------
Run numbers are counted per modality. run-2 of a T2w and run-2 of a task have
nothing to do with each other: when a subject climbs out of the scanner and
comes back, the second visit restarts every counter independently, and some
scans may not be repeated at all. Pairing by run number then coregisters an
anatomy from one head position onto functional data from another.

What does say which scans belong together is when they were acquired. This
groups every scan in a session into blocks separated by more than
--gap-minutes, picks the block the functional runs live in, and suggests the
run of each other scan from that same block. The acquisition times are written
into the file too, so the suggestion can be checked rather than trusted.
"""

import argparse
import json
import os
import re
import sys

NIFTI_EXTS = ('.nii', '.nii.gz')


# --------------------------------------------------------------------------
#  BIDS filenames
# --------------------------------------------------------------------------

def strip_ext(name):
    for ext in ('.nii.gz', '.nii', '.json'):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def entities(name):
    """{'sub': ..., 'ses': ..., 'task': ..., 'run': ...} plus 'suffix'."""
    stem = strip_ext(name)
    out = {'suffix': ''}
    for part in stem.split('_'):
        m = re.fullmatch(r'([A-Za-z0-9]+)-(.+)', part)
        if m:
            out[m.group(1).lower()] = m.group(2)
        else:
            out['suffix'] = part
    return out


def run_label(name):
    ent = entities(name)
    return 'run-%s' % ent['run'] if 'run' in ent else 'no-run'


def sidecar(path):
    base = strip_ext(path)
    js = base + '.json'
    return js if os.path.isfile(js) else None


def acq_seconds(path):
    """Seconds since midnight from the JSON sidecar, or None."""
    js = sidecar(path)
    if not js:
        return None
    try:
        with open(js) as fh:
            d = json.load(fh)
    except Exception:
        return None

    raw = d.get('AcquisitionDateTime') or d.get('AcquisitionTime')
    if not raw:
        return None
    raw = str(raw).strip()

    m = re.search(r'(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)', raw)
    if not m:
        m = re.fullmatch(r'(\d{2})(\d{2})(\d{2}(?:\.\d+)?)', raw)
    if not m:
        return None
    h, mi, se = m.groups()
    return int(h) * 3600 + int(mi) * 60 + float(se)


def series_number(path):
    js = sidecar(path)
    if not js:
        return None
    try:
        with open(js) as fh:
            return json.load(fh).get('SeriesNumber')
    except Exception:
        return None


def clock(secs):
    if secs is None:
        return '--:--:--'
    return '%02d:%02d:%02d' % (secs // 3600, (secs % 3600) // 60, secs % 60)


# --------------------------------------------------------------------------
#  Collecting one session's scans
# --------------------------------------------------------------------------

class Scan(object):
    def __init__(self, path, kind):
        self.path = path
        self.name = os.path.basename(path)
        self.kind = kind                      # the 'type' column
        self.run = run_label(self.name)
        self.acq = acq_seconds(path)
        self.series = series_number(path)
        self.block = None

    def __repr__(self):
        return '<Scan %s %s %s>' % (self.kind, self.run, clock(self.acq))


def list_nifti(folder):
    if not os.path.isdir(folder):
        return []
    out = []
    seen = set()
    for name in sorted(os.listdir(folder)):
        if name.startswith('._'):
            continue
        if not name.endswith(NIFTI_EXTS):
            continue
        stem = strip_ext(name)
        if stem in seen:          # foo.nii wins over foo.nii.gz
            continue
        seen.add(stem)
        out.append(os.path.join(folder, name))
    return out


def collect_session(ses_dir, sub, ses, reverse_pattern, fieldmap_pattern,
                    magnitude_pattern):
    """Every scan in one session, tagged with the type it will be selected as."""
    scans = []

    anat = os.path.join(ses_dir, 'anat')
    func = os.path.join(ses_dir, 'func')
    fmap = os.path.join(ses_dir, 'fmap')

    for path in list_nifti(anat):
        ent = entities(os.path.basename(path))
        if ent.get('sub') != sub[4:] or ent.get('ses') != ses[4:]:
            continue
        if ent['suffix'] in ('T1w', 'T2w'):
            scans.append(Scan(path, ent['suffix']))

    for path in list_nifti(func):
        name = os.path.basename(path)
        ent = entities(name)
        if ent.get('sub') != sub[4:] or ent.get('ses') != ses[4:]:
            continue
        if ent['suffix'] != 'bold':
            continue
        # The reverse-PE EPI is selected on its own, not as a task
        if reverse_pattern and reverse_pattern in name:
            scans.append(Scan(path, 'reverse'))
        elif ent.get('task') == 'reverse':
            scans.append(Scan(path, 'reverse'))
        elif 'task' in ent:
            scans.append(Scan(path, 'task-%s' % ent['task']))

    for path in list_nifti(fmap):
        name = os.path.basename(path)
        ent = entities(name)
        if ent.get('sub') != sub[4:] or ent.get('ses') != ses[4:]:
            continue
        suffix = ent['suffix']
        if reverse_pattern and reverse_pattern in name:
            scans.append(Scan(path, 'reverse'))
        elif suffix == 'epi':
            scans.append(Scan(path, 'reverse'))
        elif suffix in ('phasediff', 'magnitude1', 'magnitude2', 'phase1', 'phase2'):
            scans.append(Scan(path, suffix))
        elif fieldmap_pattern and strip_ext(name).endswith(fieldmap_pattern.lstrip('_')) \
                or suffix == 'fieldmap':
            scans.append(Scan(path, 'fieldmap'))
        elif magnitude_pattern and suffix == magnitude_pattern.lstrip('_'):
            scans.append(Scan(path, 'magnitude'))

    return scans


# --------------------------------------------------------------------------
#  Grouping into visits to the scanner
# --------------------------------------------------------------------------

def assign_blocks(scans, gap_sec):
    """Label each scan with a block letter; scans with no time get None.

    Blocks are runs of scans where consecutive acquisitions are no further
    apart than gap_sec. A subject leaving the scanner and coming back leaves a
    gap far larger than any within-protocol pause.
    """
    timed = sorted([s for s in scans if s.acq is not None], key=lambda s: s.acq)
    if not timed:
        return 0

    block = 0
    timed[0].block = block
    for prev, cur in zip(timed, timed[1:]):
        if cur.acq - prev.acq > gap_sec:
            block += 1
        cur.block = block
    return block + 1


def block_name(idx):
    return '?' if idx is None else chr(ord('A') + idx)


def reference_block(scans):
    """The block the functional runs are in — what everything else must match.

    The task runs are the data being analysed, so they define the visit that
    matters. When they are spread over more than one block the last one wins,
    and the caller warns: nothing can match both.
    """
    task_blocks = [s.block for s in scans
                   if s.kind.startswith('task-') and s.block is not None]
    if task_blocks:
        return max(set(task_blocks), key=task_blocks.count)
    known = [s.block for s in scans if s.block is not None]
    return max(known) if known else None


# --------------------------------------------------------------------------
#  Reporting
# --------------------------------------------------------------------------

def by_kind(scans):
    order = []
    groups = {}
    for s in scans:
        if s.kind not in groups:
            groups[s.kind] = []
            order.append(s.kind)
        groups[s.kind].append(s)
    return [(k, groups[k]) for k in sorted(order)]


def suggest(runs, ref_block):
    """Which run to use, and why. '' when nothing can be suggested safely."""
    if len(runs) == 1:
        return runs[0].run, 'only one run'

    if ref_block is not None:
        in_block = [r for r in runs if r.block == ref_block]
        if len(in_block) == 1:
            return in_block[0].run, 'the only one acquired with the functional runs'
        if len(in_block) > 1:
            # Several in the right block: the last one is the usual convention
            # (a repeat replaces the scan it repeats), but say that it is a guess.
            return in_block[-1].run, 'last of %d in the same block as the functional runs' % len(in_block)
        return '', 'NONE was acquired with the functional runs — choose by hand'

    return '', 'no acquisition times — choose by hand'


def detail(runs):
    bits = []
    for r in runs:
        stamp = clock(r.acq)
        series = '' if r.series is None else '#%d' % int(r.series)
        bits.append('%s@%s[%s]%s' % (r.run, stamp, block_name(r.block), series))
    return ','.join(bits)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bids-root', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--gap-minutes', type=float, default=20.0)
    ap.add_argument('--reverse-pattern', default='')
    ap.add_argument('--fieldmap-pattern', default='_fieldmap')
    ap.add_argument('--magnitude-pattern', default='_magnitude')
    ap.add_argument('--list', default='', help='subses_list.txt to restrict the scan to')
    ap.add_argument('--keep-existing', action='store_true',
                    help='keep the selections already in the output file')
    args = ap.parse_args()

    gap_sec = args.gap_minutes * 60

    # Selections already made, so a re-scan does not throw away your decisions
    previous = {}
    if args.keep_existing and os.path.isfile(args.output):
        previous = read_previous(args.output)

    pairs = []
    if args.list and os.path.isfile(args.list):
        with open(args.list) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    pairs.append((parts[0], parts[1]))
    else:
        for sub in sorted(os.listdir(args.bids_root)):
            if not sub.startswith('sub-'):
                continue
            sub_dir = os.path.join(args.bids_root, sub)
            if not os.path.isdir(sub_dir):
                continue
            for ses in sorted(os.listdir(sub_dir)):
                if ses.startswith('ses-') and os.path.isdir(os.path.join(sub_dir, ses)):
                    pairs.append((sub, ses))

    rows = []
    counts = {}
    split_sessions = []
    unmatched = []

    for sub, ses in pairs:
        ses_dir = os.path.join(args.bids_root, sub, ses)
        if not os.path.isdir(ses_dir):
            continue

        scans = collect_session(ses_dir, sub, ses, args.reverse_pattern,
                                args.fieldmap_pattern, args.magnitude_pattern)
        if not scans:
            continue

        n_blocks = assign_blocks(scans, gap_sec)
        ref = reference_block(scans)

        groups = by_kind(scans)
        multi = [(kind, runs) for kind, runs in groups if len(runs) > 1]
        if not multi:
            continue

        if n_blocks > 1:
            split_sessions.append((sub, ses, n_blocks))
            print_session(sub, ses, groups, ref, n_blocks)

        for kind, runs in multi:
            runs = sorted(runs, key=lambda r: (r.acq if r.acq is not None else 1e12, r.run))
            sel, why = suggest(runs, ref)

            key = (sub, ses, kind)
            if key in previous and previous[key]:
                sel = previous[key]
                why = 'kept from the previous run_selection.tsv'
            if not sel:
                unmatched.append((sub, ses, kind))

            rows.append([sub, ses, kind,
                         ','.join(r.run for r in runs),
                         sel,
                         block_name(ref),
                         detail(runs),
                         why])
            counts[kind] = counts.get(kind, 0) + 1

    write_tsv(args.output, rows)
    summarise(rows, counts, split_sessions, unmatched, args)
    return 0


def print_session(sub, ses, groups, ref, n_blocks):
    print('  %s / %s — %d separate blocks of acquisition:' % (sub, ses, n_blocks))
    flat = []
    for _kind, runs in groups:
        flat.extend(runs)
    for s in sorted(flat, key=lambda s: (s.acq if s.acq is not None else 1e12, s.name)):
        marker = ' <- functional block' if s.block == ref and s.kind.startswith('task-') else ''
        print('      [%s] %s  %-14s %s%s' % (block_name(s.block), clock(s.acq),
                                             s.kind, s.run, marker))
    print('')


def read_previous(path):
    out = {}
    with open(path) as fh:
        header = None
        for i, raw in enumerate(fh):
            row = raw.rstrip('\n').split('\t')
            if i == 0 and row and row[0].strip().lower() in ('subject', 'sub'):
                header = [c.strip().lower() for c in row]
                continue
            cols = header or ['subject', 'session', 'type', 'available_runs', 'selected_run']
            def get(name):
                if name not in cols:
                    return ''
                idx = cols.index(name)
                return row[idx].strip() if idx < len(row) else ''
            if get('subject') and get('selected_run'):
                out[(get('subject'), get('session'), get('type'))] = get('selected_run')
    return out


def write_tsv(path, rows):
    header = ['subject', 'session', 'type', 'available_runs', 'selected_run',
              'functional_block', 'runs_detail', 'notes']
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join(r) + '\n')


def summarise(rows, counts, split_sessions, unmatched, args):
    if not rows:
        print('No multi-run cases found. Every session has a single run of each scan.')
        print('')
        print('No run_selection.tsv is needed — remove RUN_SELECTION_FILE from the')
        print('config, or leave it pointing at a file that does not exist.')
        return

    print('Multi-run cases found:')
    for kind in sorted(counts):
        print('  %-18s %d subject-session(s)' % (kind, counts[kind]))
    print('')
    print('  Total rows:        %d' % len(rows))
    print('  Split sessions:    %d  (subject left the scanner between sequences)'
          % len(split_sessions))
    print('')
    print('Written to: %s' % args.output)
    print('')

    if unmatched:
        print('!! %d row(s) could NOT be suggested automatically:' % len(unmatched))
        for sub, ses, kind in unmatched[:12]:
            print('     %s %s %s' % (sub, ses, kind))
        if len(unmatched) > 12:
            print('     ... and %d more' % (len(unmatched) - 12))
        print('   Fill in their selected_run column by hand. With')
        print('   STRICT_RUN_MATCHING=true the preprocessing stops on these rather')
        print('   than guessing.')
        print('')

    print('WHAT TO DO NEXT')
    print('---------------')
    print('1. Open the file (a spreadsheet opens it as tab-separated).')
    print('2. Check the selected_run column. It is pre-filled with the run')
    print('   acquired together with the functional data; runs_detail shows every')
    print('   run with its acquisition time and its block letter [A], [B], ...')
    print('3. Correct anything the suggestion got wrong, and fill in the blanks.')
    print('4. Point RUN_SELECTION_FILE in pipeline_config.cfg at this file.')
    print('')
    print('Re-running this scanner keeps the choices you have already made')
    print('(--keep-existing, which run_pipeline.sh passes by default).')


if __name__ == '__main__':
    sys.exit(main())
