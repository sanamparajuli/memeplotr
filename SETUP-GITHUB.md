# Putting this on GitHub and shinyapps.io

Everything below assumes you are in this directory.

## 1. Replace the placeholder

`YOURNAME` appears in `DESCRIPTION` (URL, BugReports) and `README.md` (install
line, CI badge). Replace it, and put your own name and email in the
`Authors@R` field of `DESCRIPTION` while you are there — the placeholder is
the remaining `R CMD check` NOTE.

```bash
# macOS
sed -i '' 's/YOURNAME/your-github-handle/g' DESCRIPTION README.md
# Linux
sed -i    's/YOURNAME/your-github-handle/g' DESCRIPTION README.md
```

## 2. Create the repo and push

Create an empty **public** repository called `memeplotr` on GitHub — no
README, no .gitignore, no licence, since this directory already has them.
Then:

```bash
git add -A
git commit -m "memeplotr 0.1.4"
git branch -M main
git remote add origin https://github.com/your-github-handle/memeplotr.git
git push -u origin main
```

The R-CMD-check workflow runs on that push; the badge in `README.md` turns
green when it passes.

## 3. Install from GitHub

This step matters. Installing from the tarball produces a package the
deployment server cannot find; installing from GitHub records where it came
from.

```r
remotes::install_github("your-github-handle/memeplotr", build_vignettes = TRUE)
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
remotes::install_github("your-github-handle/memeplotr")
rsconnect::deployApp("deploy", appName = "memeplotr")
```

## If you would rather not use GitHub

`memeplotr::bundle_shiny_app("memeplotr-shinyapp")` writes a directory that
carries the package source inside it and loads it with `pkgload::load_all()`,
so the server only installs CRAN packages. Deploy that directory instead.
