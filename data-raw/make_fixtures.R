## ---------------------------------------------------------------------------
## Builds the simulated example files shipped in inst/extdata.
##
## These are NOT the output of a real MEME run: the sequences, motifs, PWMs,
## p-values and E-values are simulated with the generator below so that the
## package has a small, redistributable, fully reproducible example. The XML
## structure follows the MEME Suite schema (a 5.x file and a 4.x-layout file)
## so that it exercises the real parser paths.
##
##   Rscript data-raw/make_fixtures.R
## ---------------------------------------------------------------------------
set.seed(20260115)

outdir <- file.path("inst", "extdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
AA_COL <- c(A="CCFF00", C="FFFF00", D="FF0000", E="FF0066", F="00FF66", G="FF9900",
            H="0066FF", I="66FF00", K="6600FF", L="33FF00", M="00FF00", N="CC00FF",
            P="FFCC00", Q="FF00CC", R="0000FF", S="FF3300", T="FF6600", V="99FF00",
            W="00CCFF", Y="00FFCC")
BG <- c(A=.074,C=.025,D=.054,E=.054,F=.047,G=.074,H=.026,I=.068,K=.058,L=.099,
        M=.025,N=.045,P=.039,Q=.034,R=.052,S=.057,T=.051,V=.073,W=.013,Y=.032)
BG <- BG / sum(BG)

## ---- sequences and phylogeny ----------------------------------------------
tips <- c("Ath_MYB1","Ath_MYB2","Osa_MYB1","Osa_MYB2","Zma_MYB1","Zma_MYB3",
          "Sbi_MYB2","Vvi_MYB4","Ptr_MYB4","Gma_MYB5","Sly_MYB6","Stu_MYB6")
newick <- paste0(
  "((((Ath_MYB1:0.052,Ath_MYB2:0.061):0.084,(Osa_MYB1:0.073,Osa_MYB2:0.049):0.091):0.118,",
  "((Zma_MYB1:0.058,Zma_MYB3:0.041):0.072,Sbi_MYB2:0.094):0.106):0.151,",
  "(((Vvi_MYB4:0.081,Ptr_MYB4:0.069):0.063,Gma_MYB5:0.102):0.088,",
  "(Sly_MYB6:0.121,Stu_MYB6:0.113):0.079):0.134);")
writeLines(newick, file.path(outdir, "example_tree.nwk"))

seq_len_v <- c(Ath_MYB1=312, Ath_MYB2=298, Osa_MYB1=341, Osa_MYB2=327, Zma_MYB1=355,
               Zma_MYB3=288, Sbi_MYB2=364, Vvi_MYB4=241, Ptr_MYB4=255, Gma_MYB5=276,
               Sly_MYB6=398, Stu_MYB6=385)

## ---- motifs ----------------------------------------------------------------
mot <- list(
  list(w = 29, cons = "WTPEEDELLRRLVEKHGARNWSLIARSIP"),
  list(w = 21, cons = "GKSCRLRWLNYLRPDLKRGNI"),
  list(w = 15, cons = "SEEEDRIIRLHSLLG"),
  list(w = 18, cons = "NYWNTHLKKRLLSQGIDP"),
  list(w = 25, cons = "MSDLLQAPHCFSPTTEEAWMADFLN"),
  list(w = 12, cons = "KQLERSSILSSY")
)
nmot <- length(mot)

make_pwm <- function(cons) {
  w <- nchar(cons); cs <- strsplit(cons, "")[[1]]
  t(vapply(seq_len(w), function(i) {
    p <- stats::rgamma(20, shape = 0.12); names(p) <- AA
    p[cs[i]] <- p[cs[i]] + stats::runif(1, 2.0, 9.0)
    p <- p / sum(p)
    p <- round(p, 6); p[cs[i]] <- p[cs[i]] + (1 - sum(p))
    p
  }, numeric(20)))
}
pwms <- lapply(mot, function(m) make_pwm(m$cons))

regexp_of <- function(p) {
  paste(vapply(seq_len(nrow(p)), function(i) {
    o <- sort(p[i, ], decreasing = TRUE)
    k <- max(1, sum(cumsum(o) < 0.8) + 1)
    if (k == 1) names(o)[1] else paste0("[", paste(names(o)[seq_len(k)], collapse = ""), "]")
  }, character(1)), collapse = "")
}

## ---- which motifs occur in which sequence (clade-structured) ---------------
cladeA <- tips[1:7]; cladeB <- tips[8:12]
present <- matrix(FALSE, length(tips), nmot, dimnames = list(tips, paste0("motif_", 1:nmot)))
present[, 1] <- TRUE
present[, 2] <- TRUE
present[cladeA, 3] <- TRUE
present[c("Ath_MYB1","Osa_MYB1","Zma_MYB1","Sbi_MYB2"), 4] <- TRUE
present[cladeB, 5] <- TRUE
present[c("Vvi_MYB4","Ptr_MYB4","Sly_MYB6","Stu_MYB6"), 6] <- TRUE
present["Gma_MYB5", 3] <- TRUE          # a single cross-clade retention
present["Zma_MYB3", 2] <- FALSE         # one loss

## ---- site coordinates (0-based, as MEME writes them) ----------------------
anchor <- c(0.10, 0.22, 0.42, 0.56, 0.05, 0.78)
sites <- do.call(rbind, lapply(tips, function(s) {
  L <- seq_len_v[[s]]
  idx <- which(present[s, ])
  rows <- lapply(idx, function(j) {
    w <- mot[[j]]$w
    pos <- round(anchor[j] * L + stats::runif(1, -8, 8))
    pos <- max(0, min(pos, L - w))
    data.frame(sequence = s, motif = j, position = pos,
               pvalue = signif(10^stats::runif(1, -28, -9), 3),
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d[order(d$position), ]
}))
## a couple of tandem duplications, as real families show
sites <- rbind(sites,
  data.frame(sequence = "Sly_MYB6", motif = 6, position = 330, pvalue = 3.1e-11),
  data.frame(sequence = "Sbi_MYB2", motif = 1, position = 250, pvalue = 8.4e-16))
sites <- sites[order(match(sites$sequence, tips), sites$position), ]
rownames(sites) <- NULL

seq_pv <- signif(10^stats::runif(length(tips), -95, -40), 3); names(seq_pv) <- tips

## ---- XML emitters -----------------------------------------------------------
esc <- function(x) gsub("&", "&amp;", x, fixed = TRUE)
fmt <- function(v) formatC(v, format = "f", digits = 6)

array_xml <- function(p, prefix, indent) {
  paste0(indent, "<alphabet_array>\n",
         paste0(indent, "  <value letter_id=\"", prefix, AA, "\">", p[AA], "</value>",
                collapse = "\n"),
         "\n", indent, "</alphabet_array>")
}

build_xml <- function(version = "5.5.5") {
  v4 <- substr(version, 1, 1) == "4"
  pre <- if (v4) "letter_" else ""
  o <- c()
  o <- c(o, "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>")
  o <- c(o, sprintf("<MEME version=\"%s\" release=\"Simulated example for the memeplotr R package\">", version))
  o <- c(o, if (v4)
    "<training_set datafile=\"myb_family.faa\" length=\"12\">"
    else "<training_set primary_sequences=\"myb_family.faa\" primary_count=\"12\" primary_positions=\"3840\">")
  o <- c(o, if (v4)
    "  <alphabet id=\"protein\" length=\"20\">"
    else "  <alphabet name=\"Protein\" like=\"protein\">")
  o <- c(o, paste0("    <letter id=\"", pre, AA, "\" symbol=\"", AA,
                   "\" colour=\"", AA_COL[AA], "\"/>"))
  o <- c(o, "  </alphabet>")
  o <- c(o, sprintf("  <sequence id=\"sequence_%d\" name=\"%s\" length=\"%d\" weight=\"1.000000\"/>",
                    seq_along(tips) - 1L, tips, seq_len_v[tips]))
  o <- c(o, "  <letter_frequencies>", array_xml(BG, pre, "    "), "  </letter_frequencies>")
  o <- c(o, "</training_set>")

  o <- c(o, "<model>")
  o <- c(o, paste0("  <command_line>meme myb_family.faa -protein -oc . -nostatus -mod zoops -nmotifs 6 ",
                   "-minw 6 -maxw 30 -objfun classic -markov_order 0</command_line>"))
  o <- c(o, "  <host>simulated</host>", "  <type>zoops</type>",
         sprintf("  <nmotifs>%d</nmotifs>", nmot),
         "  <evalue_threshold>inf</evalue_threshold>",
         "  <object_function>E-value of product of p-values</object_function>",
         "  <minsites>2</minsites>", "  <maxsites>12</maxsites>",
         "  <strands>none</strands>",
         "  <background_frequencies source=\"--sequences--\">",
         array_xml(BG, pre, "    "), "  </background_frequencies>")
  o <- c(o, "</model>")

  o <- c(o, "<motifs>")
  for (j in seq_len(nmot)) {
    p <- pwms[[j]]; w <- mot[[j]]$w
    ev <- signif(10^(-62 + 7 * (j - 1)), 3)
    attrs <- if (v4) {
      sprintf(paste0("id=\"motif_%d\" name=\"%d\" width=\"%d\" sites=\"%d\" ic=\"%s\" re=\"%s\"",
                     " llr=\"%d\" e_value=\"%s\" bayes_threshold=\"%s\" elapsed_time=\"%s\""),
              j, j, w, sum(present[, j]), fmt(sum(.pwm_ic_local(p))), fmt(sum(.pwm_ic_local(p)) / w),
              as.integer(round(42 * w)), format(ev, scientific = TRUE), fmt(28 + j), fmt(0.3 * j))
    } else {
      sprintf(paste0("id=\"motif_%d\" name=\"%s\" alt=\"MEME-%d\" width=\"%d\" sites=\"%d\"",
                     " ic=\"%s\" re=\"%s\" llr=\"%d\" p_value=\"%s\" e_value=\"%s\"",
                     " bayes_threshold=\"%s\" elapsed_time=\"%s\""),
              j, mot[[j]]$cons, j, w, sum(present[, j]),
              fmt(sum(.pwm_ic_local(p))), fmt(sum(.pwm_ic_local(p)) / w),
              as.integer(round(42 * w)), format(ev, scientific = TRUE),
              format(ev, scientific = TRUE), fmt(28 + j), fmt(0.3 * j))
    }
    o <- c(o, sprintf("<motif %s>", attrs))
    o <- c(o, "  <scores>", "    <alphabet_matrix>")
    for (i in seq_len(w)) {
      sc <- round(100 * log2(pmax(p[i, ], 1e-5) / BG[AA]))
      o <- c(o, array_xml(stats::setNames(as.character(sc), AA), pre, "      "))
    }
    o <- c(o, "    </alphabet_matrix>", "  </scores>")
    o <- c(o, "  <probabilities>", "    <alphabet_matrix>")
    for (i in seq_len(w)) {
      o <- c(o, array_xml(stats::setNames(fmt(p[i, AA]), AA), pre, "      "))
    }
    o <- c(o, "    </alphabet_matrix>", "  </probabilities>")
    o <- c(o, paste0("  <regular_expression>\n", regexp_of(p), "\n  </regular_expression>"))
    o <- c(o, "  <contributing_sites>")
    sj <- sites[sites$motif == j, ]
    sj <- sj[!duplicated(sj$sequence), ]
    for (k in seq_len(nrow(sj))) {
      si <- match(sj$sequence[k], tips) - 1L
      cs <- strsplit(mot[[j]]$cons, "")[[1]]
      o <- c(o, sprintf("    <contributing_site sequence_id=\"sequence_%d\" position=\"%d\" strand=\"none\" pvalue=\"%s\">",
                        si, sj$position[k], format(sj$pvalue[k], scientific = TRUE)))
      o <- c(o, "      <left_flank>GSA</left_flank>", "      <site>")
      o <- c(o, paste0("        <letter_ref letter_id=\"", pre, cs, "\"/>"))
      o <- c(o, "      </site>", "      <right_flank>TPK</right_flank>",
             "    </contributing_site>")
    }
    o <- c(o, "  </contributing_sites>", "</motif>")
  }
  o <- c(o, "</motifs>")

  o <- c(o, "<scanned_sites_summary p_thresh=\"0.0001\">")
  for (s in tips) {
    si <- match(s, tips) - 1L
    ss <- sites[sites$sequence == s, ]
    o <- c(o, sprintf("  <scanned_sites sequence_id=\"sequence_%d\" pvalue=\"%s\" num_sites=\"%d\">",
                      si, format(seq_pv[[s]], scientific = TRUE), nrow(ss)))
    o <- c(o, sprintf("    <scanned_site motif_id=\"motif_%d\" strand=\"none\" position=\"%d\" pvalue=\"%s\"/>",
                      ss$motif, ss$position, format(ss$pvalue, scientific = TRUE)))
    o <- c(o, "  </scanned_sites>")
  }
  o <- c(o, "</scanned_sites_summary>", "</MEME>")
  o
}

## local copy so the script runs standalone, without loading the package
.pwm_ic_local <- function(p) {
  pp <- pmax(p, 1e-9); pp <- pp / rowSums(pp)
  log2(ncol(pp)) - (-rowSums(pp * log2(pp)))
}

writeLines(build_xml("5.5.5"), file.path(outdir, "example_meme.xml"))
writeLines(build_xml("4.12.0"), file.path(outdir, "example_meme_v4.xml"))

## ---- a small FIMO-style scan table -----------------------------------------
fimo <- do.call(rbind, lapply(seq_len(nrow(sites)), function(i) {
  data.frame(
    motif_id = paste0("motif_", sites$motif[i]),
    motif_alt_id = paste0("MEME-", sites$motif[i]),
    sequence_name = sites$sequence[i],
    start = sites$position[i] + 1L,
    stop = sites$position[i] + mot[[sites$motif[i]]]$w,
    strand = "+",
    score = round(stats::runif(1, 12, 48), 2),
    `p-value` = sites$pvalue[i],
    `q-value` = signif(sites$pvalue[i] * 25, 3),
    matched_sequence = mot[[sites$motif[i]]]$cons,
    check.names = FALSE, stringsAsFactors = FALSE)
}))
utils::write.table(fimo, file.path(outdir, "example_fimo.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)

message("wrote: ", paste(list.files(outdir), collapse = ", "))
