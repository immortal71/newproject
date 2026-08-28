import sys

COMMANDS = {"detect-circular", "score-stability", "filter-homology"}


def main() -> int:
    args = sys.argv[1:]
    if args and args[0] in COMMANDS:
        command, rest = args[0], args[1:]
    else:
        command, rest = "detect-circular", args  # backward-compatible default

    if command == "detect-circular":
        from .cli import run
    elif command == "score-stability":
        from .stability_cli import run
    else:
        from .homology_cli import run
    return run(rest)


if __name__ == "__main__":
    raise SystemExit(main())
