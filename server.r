library(dplyr)
library(DBI)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(shinyjs)

con <- dbConnect(RMySQL::MySQL(), user = 'ucsc',
         password='ucsc', host='localhost',db='ucsc_hg38')
rv <- reactiveValues()

snp_inject <- function(SEQ, CHR, START, WIDTH) {
  # default to plut strand
  # START 1 based
  iupac_code <- names(IUPAC_CODE_MAP)
  names(iupac_code) <- IUPAC_CODE_MAP
  
  query2 <- paste("SELECT chromEnd, observed, strand, alleles from snp", rv$DBVER, "Common where class = 'single'
    and chromEnd >= ", START, " and chromEnd <= ", (START + WIDTH - 1),
                  " and chrom = '", CHR, "'", sep="")
  res2 <- dbFetch(dbSendQuery(con, query2))
  code <- character()
  if (nrow(res2) > 0) {
    for (i in 1:nrow(res2)) {
      if (res2$alleles[i] != "") {
        observed <- res2$alleles[i]
      } else {
        observed <- res2$observed[i]
      }
      nt <- strsplit(x = observed, split = "\\W+", perl = T) %>% unlist %>% sort %>% paste(collapse="")
      code <- c(code, iupac_code[nt])
    }
  }
  return <- list()
  return$injected <- replaceLetterAt(x = SEQ, at = res2$chromEnd - START + 1, letter = code)
  return$poss <- res2$chromEnd - START + 1
  return(return)
}
# snp_inject(DNAString("ATCGATCGAT"), "chr1", 13111, 10)

subchar <- function(string, pos, char="-") { 
  # http://r.789695.n4.nabble.com/String-position-character-replacement-td4370354.html
  for(i in pos) { 
    string <- gsub(paste("^(.{", i-1, "}).", sep=""), paste("\\1", char, sep = ""), string) 
  } 
  string 
} 

subchar2lower <- function(string, pos) { 
  for(i in pos) { 
    string <- gsub(paste("^(.{", i-1, "})(.)", sep=""), "\\1\\L\\2", string, perl = T) 
  } 
  string 
}

function(input, output) {
  output$text <- renderPrint(IUPAC_CODE_MAP)
  
  observeEvent(input$submit, {
    rv$snp_tbl <- data.frame(SNP_ID = character(), Strand = character(), Class = character(),
                             Sequence = character())
    rs <- strsplit(input$rstext, "\\W+", perl = T) %>% unlist
    FLANK_LEN <- input$FLANK_LEN
    rv$DBVER <- input$dbsnp
    rv$FORMAT <- input$format
    
    withProgress(message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        for (i in 1:length(rs)) {
          name <- rs[i]
          query1 <- paste("select chromStart,chromEnd,chrom, strand, observed, alleles, class from snp", rv$DBVER, " where name = '", name, "' AND chrom not like '%\\_%'", sep = "")
          res <- dbFetch(dbSendQuery(con, query1))
          if (nrow(res) == 1) {
            seq_up <- subseq(x = Hsapiens[[res$chrom[1]]], start = res$chromStart[1] - FLANK_LEN + 1, width = FLANK_LEN)
            seq_dn <- subseq(x = Hsapiens[[res$chrom[1]]], start = res$chromEnd[1] + 1, width = FLANK_LEN)
            
            if (rv$FORMAT != "RAW") {
              # replace common SNP with IUPAC CODE
              inject_up <- snp_inject(SEQ = seq_up, CHR = res$chrom[1], START = res$chromStart[1] - FLANK_LEN + 1, WIDTH = FLANK_LEN)
              up_poss <- inject_up$poss
              inject_dn <- snp_inject(SEQ = seq_dn, CHR = res$chrom[1], START = res$chromEnd[1] + 1, WIDTH = FLANK_LEN)
              dn_poss <- inject_dn$poss
              if (rv$FORMAT == "IUPAC") {
                seq_up <- inject_up$inject
                seq_dn <- inject_dn$inject
              }
            }
            
            if (res$strand[1] == "-") {
              seq_tmp <- reverseComplement(seq_dn)
              seq_dn <- reverseComplement(seq_up)
              seq_up <- seq_tmp
              tmp_poss <- up_poss
              up_poss <- FLANK_LEN + 1 - dn_poss
              dn_poss <- FLANK_LEN + 1 - tmp_poss
            }
            if (rv$FORMAT == "Lower Case") {
              seq_up <- subchar2lower(seq_up, up_poss)
              seq_dn <- subchar2lower(seq_dn, dn_poss)
            }
            if (rv$FORMAT == "n") {
              seq_up <- subchar(seq_up, up_poss, "n")
              seq_dn <- subchar(seq_dn, dn_poss, "n")
            }
            if (rv$FORMAT == "IUPAC") {
              seq_up <- gsub("([MRWSYKVHDBN])", "\\L\\1", seq_up, perl = T)
              seq_dn <- gsub("([MRWSYKVHDBN])", "\\L\\1", seq_dn, perl = T)
            }
            if (res$alleles[1] != "") {
              observed <- sub(pattern = ",$", x = res$alleles[1], replacement = "")
              observed <- gsub(pattern = ",", x = observed, replacement = "/")
            } else {
              observed <- res$observed[1]
            }
            seq <- paste(seq_up, "[", observed, "]", seq_dn, sep="")
            rv$snp_tbl <- rbind(rv$snp_tbl, data.frame(SNP_ID = name, Strand = res$strand[1],
                                                 Class= res$class[1], Sequence = seq))
          }
          incProgress(1/length(rs), name)
        }
      })
    
    output$tbl <- renderTable({rv$snp_tbl})
  })
  
  observeEvent(input$highlight, {
    runjs("$('.shiny-table tr td:nth-child(4)').highlight('[atcgnmrwsykvhdb]');")
  })
    output$downloadData <- downloadHandler(
    filename = function() { 
      paste(Sys.Date(), '.txt', sep='') 
    },
    content = function(file) {
      write.table(rv$snp_tbl[,c(1,4)], file, quote = F, row.names = F, col.names = T, sep = "\t")
    }
  )
}
