options(stringsAsFactors = FALSE)

if( system("which wkhtmltopdf > /dev/null")!=0 ) stop(
         "wkhtmltopdf not installed, stop")

args <- commandArgs(trailingOnly=TRUE)
htmlfile <- args[1]

if( is.na(htmlfile) | !file.exists(htmlfile) ) stop(
    "Usage:\n Rscript printReportToPdf.R file.html")

output <- sub("html$", "pdf.html", htmlfile)

text <- readLines(htmlfile)
##text <- scan(file=htmlfile, what = "character", sep ="\n",  allowEscapes = TRUE)

i.body <- grep("<\\s*body\\s*>", text)
i.toc0 <- grep("id=.*toc_0\"", text) ## title
i.toc1 <- grep("id=.*toc_1\"", text) ## introduction

stopifnot(length(i.body)==1)
stopifnot(length(i.toc0)==1)
stopifnot(length(i.toc1)==1)
stopifnot(i.body<i.toc0)
stopifnot(i.toc0<i.toc1)

i.hr1 <- grep("<hr", text[i.toc0:(i.toc1-1)]) ## page brake
if( length(i.hr1)==1 ) i.toc1 <- i.toc0 + i.hr1 -1

text[i.body] <- paste(text[c(i.body, i.toc0:(i.toc1-1) )], collapse = "\n")
text[i.toc0:(i.toc1-1)] <- ""

i.toc.start <- grep("#toc\\s*\\{", text)
stopifnot(length(i.toc.start)==1)
i.toc.end <- i.toc.start + grep("\\}", text[i.toc.start+1:50])[1]
toccss <- text[i.toc.start:i.toc.end]

i.position <- grep("position", toccss, ignore.case=TRUE) + i.toc.start -1
i.width <- grep("width", toccss, ignore.case=TRUE) + i.toc.start -1

stopifnot( grepl("position", text[i.position], ignore.case = TRUE) )
stopifnot( grepl("width", text[i.width], ignore.case = TRUE) )

text[i.width] <- "  width: 90%;"
text[i.position] <- "  position: static;"

i.marginleft <- grep("margin-left:210px", text)
text[i.marginleft] <- "margin-left:50px;"

write(text, output)

cmd <- sprintf("wkhtmltopdf -L 2mm -s Letter %s %s", output, sub(".html$", "", output))
message(cmd)
stopifnot( system(cmd)==0 )
unlink(output)

