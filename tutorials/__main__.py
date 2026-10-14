"""Read or run the examples, starting with the graph tutorial and ABA."""
import argparse
import contextlib
import io
import json
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]


def examples():
    return json.loads((ROOT / 'tutorials/minimal_examples.json').read_text(encoding="utf-8"))


def main(argv=None):
    rows = examples()
    choices = {row['id']: row for row in rows}
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('example', nargs='?', choices=choices)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--source', nargs='?', const='tutorial',
                        choices=('tutorial', 'inputs', 'counting', 'events'),
                        help='print the actual source instead of running it')
    action.add_argument('--out', type=Path, help='save a readable summary, observable definitions and any closed ODE')
    parser.add_argument('--details', action='store_true', help='also save the complete JSON record and available traces')
    args = parser.parse_args(argv)
    if args.details and (args.out is None or args.source is not None):
        parser.error('--details requires --out and an example execution')
    if args.out is not None:
        if args.example is None:
            parser.error('--out requires an example')
        if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
            parser.error(f'{args.out} must be a new or empty output directory')
    if args.example is None:
        for row in rows:
            print(f"{row['id']:16} {row['title']}")
        print('\nRun:  python -m tutorials aba')
        print('Read: python -m tutorials aba --source')
        print('Next: python -m tutorials aba --source counting')
        return 0
    row = choices[args.example]
    if args.source:
        name = ('tutorials/' + row['script'] + '.py' if args.source == 'tutorial'
                else row['sources'].get(args.source))
        if name is None:
            parser.error(f"{args.example} has no separate {args.source} source")
        print((ROOT / name).read_text(encoding="utf-8"), end='')
    else:
        # The scripts remain usable directly; the dispatcher uses one result formatter.
        with contextlib.redirect_stdout(io.StringIO()):
            namespace = runpy.run_module('tutorials.' + row['script'], run_name='__main__')
        from mfobs.formatting import render
        result = namespace['result']
        if hasattr(result, 'record'):
            result = result.record()
        print(render(result, name=row['title']))
        if args.out is not None:
            from mfobs.outputs import export_summary
            export_summary(result, args.out, name=row['title'], details=args.details)
            print(f'\nSaved results: {args.out.resolve()} (start with summary.txt)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
