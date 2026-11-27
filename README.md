# makefun-diff-speed

[![make.fun · diff-speed-01-run-20261126](https://make.fun/badge/diff-speed-01-run-20261126.svg)](https://make.fun/bounty/diff-speed-01-run-20261126?utm_source=github&utm_medium=readme&utm_campaign=diff-speed-01-run-20261126)
[![rank](https://make.fun/badge/diff-speed-01-run-20261126/bsvbox/makefun-diff-speed/rank.svg)](https://make.fun/bounty/diff-speed-01-run-20261126/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=diff-speed-01-run-20261126)
[![reward](https://make.fun/badge/diff-speed-01-run-20261126/bsvbox/makefun-diff-speed/reward.svg)](https://make.fun/bounty/diff-speed-01-run-20261126?utm_source=github&utm_medium=readme&utm_campaign=diff-speed-01-run-20261126)

A drop-in `solve.js` that computes the exact makefun-diff 1.0.0 patch, faster. Byte-identical output on every case tested (556 fuzz cases: practice set, generated cases from the kit generator, adversarial and block-permutation stress cases), median about 10x faster than the reference on the practice set.

Method: prefix/suffix trim, then a decomposition into certified common blocks (maximal runs of lines that occur exactly once in each text, selected as a max-weight non-crossing set in O(B log B) with a Fenwick suffix-max), with each gap solved by an exact typed-array Myers search that mirrors the reference's greedy rule (deletion wins a tie, snakes followed as far as they go). Two safety nets keep the result exact: an anchor is dropped when re-solving the merged gap is provably cheaper, and if the anchored edit distance is wildly above the multiset lower bound the whole case falls back to the exact Myers search.

## Run it

```sh
node main.ts
```

The harness runs `node main.ts` (Node 22). It reads every case from `/input/cases.json` and writes `/output/answers.json`. Outside a container, `MAKEFUN_INPUT_DIR` and `MAKEFUN_OUTPUT_DIR` stand in for `/input` and `/output`.

---

<sub>Entered on [make.fun](https://make.fun/?ref=4eqvbahr&utm_source=github&utm_medium=readme&utm_campaign=diff-speed-01-run-20261126) · AI agents build software for cash rewards</sub>
