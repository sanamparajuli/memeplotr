# Putting this on GitHub and shinyapps.io

Everything below assumes you are in this directory.

## 1. Maintainer details

The GitHub handle `sanamparajuli` is already substituted into `DESCRIPTION`
(URL, BugReports) and `README.md` (install line, CI badge).

One placeholder is left — the maintainer in `Authors@R`, which is the
remaining `R CMD check` NOTE. Edit `DESCRIPTION` line 6 to your real name and
email:

```r
Authors@R: c(
    person("Sanam", "Parajuli", email = "you@real.address", role = c("aut", "cre")))
```

Add `comment = c(ORCID = "0000-0002-...")` if you have an ORCID; leave it out
otherwise, since an invalid one raises its own NOTE.

## 2. Create the repo and push

Create an empty **public** repository called `memeplotr` on GitHub — no
README, no .gitignore, no licence, since this directory already has them.
Then:

```bash
git add -A
git commit -m "memeplotr 0.1.4"
git branch -M main
git remote add origin https://github.com/sanamparajuli/memeplotr.git
git push -u origin main
```

The R-CMD-check workflow runs on that push; the badge in `README.md` turns
green when it passes.

## 3. Install from GitHub

This step matters. Installing from the tarball produces a package the
deployment server cannot find; installing from GitHub records where it came
from.

```r
remotes::install_github("sanamparajuli/memeplotr", build_vignettes = TRUE)
packageDescription("memeplotr")$RemoteType   # must print "github"
```

If that prints `NULL`, run `remove.packages("memeplotr")` and install again
from GitHub.

## 4. Deploy

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<account>", token = "<token>",
                          secret = "<secret>")   # Account > Tokens on shinyapps.io

shiny::runApp("deploy")                          # check locally first
rsconnect::deployApp("deploy", appName = "memeplotr")
```

Your app is at `https://<account>.shinyapps.io/memeplotr/`. The first deploy
takes several minutes while the server builds ggplot2, ape and friends.

## Redeploying after a change

Push first, then reinstall, then deploy — in that order. The server installs
whatever commit is on GitHub, so deploying without pushing silently ships the
previous version.

```bash
git add -A && git commit -m "..." && git push
```
```r
remotes::install_github("sanamparajuli/memeplotr")
rsconnect::deployApp("deploy", appName = "memeplotr")
```

## If you would rather not use GitHub

`memeplotr::bundle_shiny_app("memeplotr-shinyapp")` writes a directory that
carries the package source inside it and loads it with `pkgload::load_all()`,
so the server only installs CRAN packages. Deploy that directory instead.
