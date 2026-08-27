# README `--help` block: screenshot → generated code block

Reference for repeating this change in the other bixBeta pipelines (`10x`,
`atac`, `smrna`, `nextflow`, …). They share the same flat layout — `params.*`,
`if(params.help){ log.info """…""" }`, and a banner at file scope — so the same
steps apply. Point a session at this file by absolute path.

## The problem

A README that embeds `img/success-*.png` showing terminal output goes stale
silently. The screenshot pins things that change on every run and every machine:

- the Nextflow version banner and its "a new version is available" nag
- `revision: 4c7bbc5899` — the commit it happened to run at
- the hostname in the shell prompt
- **every flag that existed the day it was taken**

In `sartools` the screenshot was missing `--trimmer` entirely, and the `g2`
screenshot showed `Usage: … -r main` while the README two lines above said
`-r g2`. Nobody noticed because an image cannot be grepped or diffed.

## The fix, in order

### 1. Make the usage line self-correcting first

Do this **before** touching the README, or you bake a branch name into docs on
every branch. In `main.nf`:

```groovy
    nextflow run https://github.com/<org>/<repo> -r ${workflow.revision ?: 'main'} < args ... >
```

`workflow.revision` is the branch/tag from `-r`. It is **null** for a local
`nextflow run main.nf`, so the elvis fallback is required — without it local
runs print `-r null`. Verified: `-r main` → `main`, `-r g2` → `g2`, local →
`main`.

This keeps `main.nf` byte-identical across branches. Hardcoding `-r g2` on the
Singularity branch would work but creates a permanent one-line divergence that
every future `main` → `g2` merge has to step around.

### 2. Generate the block — never retype it

Extract from source rather than transcribing, so it cannot drift. Note the
`REV` substitution: after step 1 the source holds a Groovy expression, and the
README must show what the branch actually *prints*, not the expression.

```python
import re

REV   = "main"          # the branch being edited: "main", "g2", ...
IMAGE = "![](img/success-main.png)"

src = open("main.nf").read()
help_txt = re.search(r'^log\.info """\n(.*?)^"""', src, re.S | re.M).group(1).rstrip("\n")

# render the usage line the way this branch prints it
help_txt = re.sub(r"-r \$\{workflow\.revision[^}]*\}", f"-r {REV}", help_txt)

readme = open("README.md").read()
assert IMAGE in readme, "image line not found"
open("README.md", "w").write(readme.replace(IMAGE, "```\n" + help_txt + "\n```", 1))
```

### 3. Confirm against real output, not against the source

The block is *rendered* output, so byte-comparing it to `main.nf` will fail on
the usage line by design. Compare it to an actual `--help` run instead (see
[Verifying without Nextflow installed](#verifying-without-nextflow-installed)),
normalising the `-r` token the same way the CI step does.

### 4. Include only the banner, not the Nextflow preamble

Drop `N E X T F L O W ~ version …`, the update nag, and the `Launching … [g2]`
line. That preamble is the part that ages worst — it is machine- and
run-specific.

### 5. Delete the orphaned images

`git rm` the screenshot and any already-unreferenced siblings. Check nothing
else points at them first; in these repos the qmd templates reference a *remote*
logo (`raw.githubusercontent.com/bixBeta/atac/main/img/trex-mini.png`), which is
unaffected by deleting a local `img/`.

Note a published release tagged before this change still renders the old README
with a now-broken image link. Cut a patch tag if that matters.

## The guardrail

Add a CI step so it cannot rot again. Requires `--help` output captured to
`/tmp/help.txt` in an earlier step.

```yaml
      # The README code block is the rendered --help output. It went stale as a
      # screenshot before; this keeps it honest. The -r revision differs by
      # design (workflow.revision is null locally, so --help prints the "main"
      # fallback while the README shows the branch), so normalise that token.
      - name: README help block matches actual output
        run: |
          set -euo pipefail
          sed -n '/^S A R  T O O L S/,/--trimmer/p' /tmp/help.txt \
            | sed -E 's@-r [A-Za-z0-9._-]+ < args@-r REV < args@' > /tmp/a.txt
          sed -n '/^S A R  T O O L S/,/--trimmer/p' README.md \
            | sed -E 's@-r [A-Za-z0-9._-]+ < args@-r REV < args@' > /tmp/b.txt
          diff -u /tmp/b.txt /tmp/a.txt || {
            echo "README help block is out of date - repaste 'nextflow run . --help'"; exit 1; }
```

Adjust the two `sed -n` ranges to the repo's banner text and last flag.

**Normalising `-r` is load-bearing.** CI runs `nextflow run . --help` on a local
checkout, where `workflow.revision` is null, so the output says `-r main` while
a g2 README says `-r g2`. Without the substitution this fails on every branch
except main.

## Gotchas

| Gotcha | Consequence |
|---|---|
| `workflow.revision` is null on local runs | Prints `-r null` without the `?:` fallback |
| GitHub Actions does not support YAML anchors | `&code`/`*code` parses in a YAML loader but Actions rejects it — duplicate the path list |
| `pull_request` needs its own `paths:` | A path filter on `push` alone still runs the full matrix on doc-only PRs |
| No `nextflow_schema.json` in these repos | An unknown `--flag` is accepted **silently**; a half-wired param renders a plausible report with stale values instead of erroring |
| Nextflow ≥ 26 | Strict parser rejects the file-scope statements these pipelines use — pin CI below 26 (`25.04.1` = server, `25.10.7` = latest 25.x) |
| Templates mount from `projectDir` | qmd/README changes need **no** image rebuild; `docker-build.yml` is path-filtered to `Dockerfile` anyway |
| `nextflow run <org>/<repo> -r main` | Reads `~/.nextflow/assets/…`, not your clone, and `-r main` does not re-fetch — `nextflow pull` after every push |

## Verifying without Nextflow installed

Nextflow is not on the PATH on the Mac, but Docker is:

```bash
docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w -e NXF_HOME=/w/.nxf \
  nextflow/nextflow:25.04.1 nextflow run . --help
```

Swap `run . --help` for `config .` to check the resolved config. Run against a
`git archive <branch> | tar -x -C <tmpdir>` export to avoid leaving `work/` and
`.nextflow/` in the repo.

## Checklist for the next repo

1. `grep -rn 'img/' README.md` — find the screenshot
2. Read the image to see what it claims, and diff that against real `--help`
3. Change the usage line to `${workflow.revision ?: 'main'}`
4. Run `--help` in the container; confirm the fallback prints
5. Generate the block, substituting the branch's revision
6. `git rm` the orphaned images
7. Add the CI step, adjusting the `sed` ranges
8. Push, confirm CI green, then `nextflow pull` before the next real run
