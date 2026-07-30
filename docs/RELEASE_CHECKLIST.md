# GitHub release checklist

- [ ] Confirm the manuscript title and repository name.
- [ ] Confirm the author/maintainer metadata.
- [ ] Choose a software license with the study team or institution.
- [ ] Confirm whether the repository should be public immediately or private through peer review.
- [ ] Obtain approval for every aggregate table or figure proposed for release.
- [ ] Keep all raw and individual-level files outside the repository.
- [ ] Run `Rscript run_all.R --mode=check`.
- [ ] Reproduce the full pipeline in a clean local results directory.
- [ ] Compare final aggregate estimates and figure checksums with the locked manuscript outputs.
- [ ] Initialize Git only after the above checks pass.
- [ ] Review staged filenames and file contents before the first push.
- [ ] Create a tagged release corresponding to the accepted manuscript version.

