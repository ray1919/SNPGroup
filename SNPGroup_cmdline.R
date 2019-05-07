library(dplyr)
library(DBI)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

# Update: 2017-03-07
# ! minus strand common snp iupac code bug killed
# + download nearby snp (<= 30nt) freqs in csv table
# Update: 2017-04-24
# command line version

rv <- list()
rv$DBVER <- "150"
rv$FORMAT <- "Lower Case" # c("IUPAC", "Lower Case", "n")
rv$FLANK_LEN <- 200 # 上下游获取长度
rv$dbuser <- "ucsc"
rv$dbname <- "ucsc_hg38"
rv$dbpass <- "ucsc"
rv$dbhost <- "localhost"

snp_inject <- function(SEQ, CHR, START, WIDTH) {
  # default to plut strand
  # START 1 based
  iupac_code <- names(IUPAC_CODE_MAP)
  names(iupac_code) <- IUPAC_CODE_MAP
  con <- dbConnect(RMySQL::MySQL(), user = rv$dbuser,
                   password=rv$dbpass, host=rv$dbhost,db=rv$dbname)
  
  query2 <- paste("SELECT chromEnd, observed, strand, alleles, alleleFreqs from snp",
                  rv$DBVER, "Common where class = 'single' and chromEnd >= ", START,
                  " and chromEnd <= ", (START + WIDTH - 1), " and chrom = '", CHR, "'", sep="")
  res2 <- dbFetch(dbSendQuery(con, query2))
  code <- character()
  dbDisconnect(con)
  if (nrow(res2) > 0) {
    for (i in 1:nrow(res2)) {
      if (res2$alleles[i] != "") {
        observed <- res2$alleles[i]
      } else {
        observed <- res2$observed[i]
      }
      nt <- strsplit(x = observed, split = "\\W+", perl = T) %>% unlist %>% sort %>% paste(collapse="")
      if (res2$strand[i] == "+") {
        code_i <- iupac_code[nt]
      } else {
        code_i <- iupac_code[nt] %>% DNAString %>% reverseComplement %>% toString
      }
      code <- c(code, code_i)
    }
  }
  res2$relPos <- res2$chromEnd - START + 1
  return <- list()
  return$injected <- replaceLetterAt(x = SEQ, at = res2$chromEnd - START + 1, letter = code)
  return$poss <- res2$chromEnd - START + 1
  return$relPos <- res2
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

rs2snpgroup <- function(input) {


    rv$snp_tbl <- data.frame(SNP_ID = character(), Strand = character(), Class = character(),
                             Sequence = character())
    rv$relPos <- data.frame()
    rs <- strsplit(input, "\\W+", perl = T) %>% unlist
    FLANK_LEN <- rv$FLANK_LEN
    con <- dbConnect(RMySQL::MySQL(), user = rv$dbuser,
                     password=rv$dbpass, host=rv$dbhost,db=rv$dbname)
    

    for (i in 1:length(rs)) {
      name <- rs[i]
      print(name)
      query1 <- paste("select chromStart,chromEnd,chrom, strand, observed, alleles, class from snp", rv$DBVER, " where name = '", name, "' AND chrom not like '%\\_%'", sep = "")
      res <- dbFetch(dbSendQuery(con, query1))
      if (nrow(res) == 1) {
       seq_up <- subseq(x = Hsapiens[[res$chrom[1]]], start = res$chromStart[1] - FLANK_LEN + 1, width = FLANK_LEN)
       seq_dn <- subseq(x = Hsapiens[[res$chrom[1]]], start = res$chromEnd[1] + 1, width = FLANK_LEN)
       
       if (rv$FORMAT != "RAW") {
         # replace common SNP with IUPAC CODE
         inject_up <- snp_inject(SEQ = seq_up, CHR = res$chrom[1], START = res$chromStart[1] - FLANK_LEN + 1, WIDTH = FLANK_LEN)
         up_relPos <- inject_up$relPos
         up_poss <- inject_up$poss
         inject_dn <- snp_inject(SEQ = seq_dn, CHR = res$chrom[1], START = res$chromEnd[1] + 1, WIDTH = FLANK_LEN)
         dn_relPos <- inject_dn$relPos
         dn_poss <- inject_dn$poss
         if (rv$FORMAT == "IUPAC") {
           seq_up <- inject_up$inject
           seq_dn <- inject_dn$inject
         }
       }
       if (nrow(up_relPos) > 0) {
         up_relPos$rs <- name
       }
       if (nrow(dn_relPos) > 0) {
         dn_relPos$rs <- name
       }
       if (res$strand[1] == "+") {
         up_relPos$relPos <- -(FLANK_LEN - up_relPos$relPos + 1)
       } else if (res$strand[1] == "-") {
         up_relPos$relPos <- FLANK_LEN - up_relPos$relPos + 1
         dn_relPos$relPos <- -dn_relPos$relPos
         seq_tmp <- reverseComplement(seq_dn)
         seq_dn <- reverseComplement(seq_up)
         seq_up <- seq_tmp
         tmp_poss <- up_poss
         up_poss <- FLANK_LEN + 1 - dn_poss
         dn_poss <- FLANK_LEN + 1 - tmp_poss
       }
       rv$relPos <- rbind(rv$relPos, rbind(up_relPos,dn_relPos))
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
    }
    
    dbDisconnect(con)
    
    output <- rv$snp_tbl

    seqfile <- paste("seq", Sys.Date(), '.txt', sep='')
    rv$snp_tbl$SNP_ID <- paste(rv$snp_tbl$SNP_ID, rv$snp_tbl$Strand, sep="")
    write.table(rv$snp_tbl[,c(1,4)], seqfile, quote = F, row.names = F, col.names = T, sep = "\t")

    # 30bp 以内的常见snp
    nbsnpfile <- paste("nbsnp", Sys.Date(), '.csv', sep='')
    write.csv(filter(rv$relPos, abs(relPos) < 31), nbsnpfile, quote = T, row.names = F, col.names = T, sep = ",")

    return(output)
}

table_return <- rs2snpgroup("rs1042522
rs1143627
rs16944
rs1800795
rs1800872
rs1800896
rs20417
rs20432
rs2227983
rs4073
rs4646903
rs861539")



