#!/bin/sh
# fetch.sh - download closing prices from Stooq and merge them into the
# CSV read by cvarport: header `date,T1,...,Tn`, only dates common to
# every ticker, ascending.
#
# usage: scripts/fetch.sh [-i d|w|m] TICKER... > prices.csv
#
# Tickers use Stooq notation (e.g. spy.us, eunl.de): check each one on
# stooq.com first. STOOQ_URL overrides the download endpoint.

set -eu
LC_ALL=C
export LC_ALL

usage() {
    echo "usage: $0 [-i d|w|m] TICKER..." >&2
    exit 2
}

interval=d
while getopts i: opt; do
    case $opt in
        i) interval=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))
[ $# -ge 1 ] || usage
case $interval in d|w|m) ;; *) usage ;; esac

base=${STOOQ_URL:-https://stooq.com/q/d/l/}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
trap 'exit 1' INT TERM

header=date
for t in "$@"; do
    curl -fsS "$base?s=$t&i=$interval" > "$tmp/raw.csv"
    # Stooq can answer 200 with an HTML page (a browser check, for
    # instance): accept only its CSV, Date,Open,High,Low,Close,Volume.
    if ! head -n 1 "$tmp/raw.csv" | grep -q '^Date,Open,High,Low,Close'; then
        echo "$0: $t: Stooq did not send CSV data, download it from a browser" >&2
        exit 1
    fi
    awk -F, 'NR > 1 && $5 != "" { print $1 "," $5 }' "$tmp/raw.csv" |
        sort -t, -k1,1 > "$tmp/one.csv"
    if [ ! -s "$tmp/one.csv" ]; then
        echo "$0: no data for $t" >&2
        exit 1
    fi
    if [ -f "$tmp/all.csv" ]; then
        join -t, "$tmp/all.csv" "$tmp/one.csv" > "$tmp/next.csv"
        mv "$tmp/next.csv" "$tmp/all.csv"
    else
        mv "$tmp/one.csv" "$tmp/all.csv"
    fi
    header=$header,$t
done

echo "$header"
cat "$tmp/all.csv"
