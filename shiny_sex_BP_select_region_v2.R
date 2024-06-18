library(shiny)
library(data.table)  # for fread, which is generally faster and can handle larger files
library(dplyr)

readData <- function() {

sbp1 = fread("UKB_SBP_bolt_imputed_result_MAF0.01.tsv",drop = c("A1FREQ","BETA","SE","P_BOLT_LMM","P_BOLT_LMM.GC"))
dbp1 = fread("UKB_DBP_bolt_imputed_result_MAF0.01.tsv",drop = c("A1FREQ","BETA","SE","P_BOLT_LMM","P_BOLT_LMM.GC"))
pp1 = fread("UKB_PP_bolt_imputed_result_MAF0.01.tsv",drop = c("A1FREQ","BETA","SE","P_BOLT_LMM","P_BOLT_LMM.GC"))
refgene=fread('refGene_hg19_v2.txt',header=T)
setDF(sbp1)
setDF(dbp1)
setDF(pp1)
setDF(refgene)
datasets=list("SBP"=sbp1,"DBP"=dbp1,"PP"=pp1,"ref"=refgene)
return(datasets)
}



# Create a global variable to store the data
datasets <- readData()

ui <- fluidPage(
  titlePanel("Subset Data by CHR and BP Range"),
  sidebarLayout(
    sidebarPanel(
	radioButtons("inputMethod", "Choose Input Method:",
                   choices = list("By SNP" = "SNP", "By CHR and BP Range" = "CHR_BP"),
                   selected = "CHR_BP"),
	# Conditional UI for input fields
      uiOutput("variableInputs"),
	  
      selectInput("select", h3("BP trait"), choices = list("SBP" = 'SBP', "DBP" = 'DBP',"PP" = 'PP'), selected = 'PP'),
	  actionButton("submit", "Submit", class = "btn-primary"),
	  
	# Add citation here
      tags$hr(),  # Horizontal line for separation
      tags$p("Cite our paper here:",
             tags$i("Yang, ML., Xu, C., Gupte, T. et al. Sex-specific genetic architecture of blood pressure. Nat Med 30, 818–828 (2024)."),
             style = "font-size: 12px; text-align: left; margin-top: 20px;")  
	  
    ),
    mainPanel(
      plotOutput("plot", click = "PLOT", height = "800px"),
	  downloadButton("downloadData", "Download Table"),
      div(style = "width: 100%; text-align: center;",
          h3("Sex-stratified Summary statistics", style = "font-size: 18px;"),
		  uiOutput("displayMessage"), 
          div(style = "display: inline-block; margin-left: auto; margin-right: auto;",
              tableOutput("subsetTable")
          )  # Closing div for the table
      )
	

    )
  )
)


server <- function(input, output) {

  # Dynamically render the input UI elements based on the selected method
  output$variableInputs <- renderUI({
    if (input$inputMethod == "CHR_BP") {
      tagList(
        numericInput("chr", "Chromosome (CHR):", value = 1, min = 1, max = 23),
        numericInput("bpStart", "BP Start:", value = 110546007),
        numericInput("bpEnd", "BP End:", value = 111046007)
      )
    } else if (input$inputMethod == "SNP") {
	 tagList(
      textInput("snp", "Enter SNP ID:", value = "rs123"),
	  numericInput("window","SNP +/- bp window", value=250000)
	  )
    }
  })

  # Reactive expression for the subset data
  subset_data <- reactive({
    if (input$inputMethod == "CHR_BP") {
      # Ensure that the input values are numeric and within the expected range
      chr <- as.numeric(input$chr)
      bpStart <- as.numeric(input$bpStart)
      bpEnd <- as.numeric(input$bpEnd)

      # Filter the data based on the input values
      selected_data <- datasets[[input$select]]
      selected_data %>%
        filter(CHR == chr, BP >= bpStart, BP <= bpEnd)
    } else {
      # Assume the datasets list has a way to fetch data by SNP
	  wd <- as.numeric(input$window)
      selected_data <- datasets[[input$select]]
      snp_data<-selected_data %>%
        filter(SNP == paste(input$snp)) %>%
		select(CHR,BP)
	 
	   # Check if snp_data is empty to prevent errors in the next step
	if (nrow(snp_data) == 0) {
	    showNotification("No SNP found", type = "error")
		return(data.frame())  # Return an empty data frame if no SNP matches
	}
	  chr <- snp_data$CHR[1]
	  bp <- snp_data$BP[1]
	 # Second subset based on CHR and BP window
     selected_data %>%
     filter(CHR == chr, BP >= (bp - wd), BP <= (bp + wd))
  }})
  
  # Data limited to 100 rows for display
  display_data <- reactive({
    tb=subset_data()
    head(subset_data(), 100)
  })
  
 
  # Output the subset data as a table
  output$subsetTable <- renderTable({
    display_data()
  })

# Message about the display limit
  output$displayMessage <- renderUI({
    HTML("<div style='color: red;'>Only the top 100 rows are displayed in the table. Click the download button to access the entire dataset. </div>")
  })
  
  
  # Server logic for download
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", input$select,"-CHR",input$chr,"_", input$bpStart,"_",input$bpEnd,".csv", sep = "")
    },
    content = function(file) {
      write.csv(subset_data(), file, row.names = FALSE)
    }
  )
  
  # Output the plot
  output$plot <- renderPlot({
    # Access the reactive subset_data within the renderPlot
    sd <- subset_data()
    
    # Check if the columns exist before plotting
    if("P_BOLT_LMM.male" %in% names(sd) && "P_BOLT_LMM.female" %in% names(sd)) {
		par(mfrow=c(4,1))
		max=max(c(-log10(sd$P_BOLT_LMM.male),-log10(sd$P_BOLT_LMM.female),10))
		plot(sd$BP, -log10(sd$P_BOLT_LMM.GC.male), main = "UKB Male-only Genetic Association",pch=19,col='blue',ylim=c(0,max),xlim=c(min(sd$BP),max(sd$BP)),xlab='POS(hg19)',ylab='-log10(P)')
		abline(h=7.3,col='dark gray')
		plot(sd$BP, -log10(sd$P_BOLT_LMM.GC.female), main = "UKB Female-only Genetic Association",pch=19,col='purple',ylim=c(0,max),xlim=c(min(sd$BP),max(sd$BP)),xlab='POS(hg19)',ylab='-log10(P)')
		abline(h=7.3,col='dark gray')
		plot(sd$BP, -log10(sd$diffBetaSEX.GCpvalue), main = "UKB Sex-dimorphic effect of genetic Association",pch=19,col='dark orange',ylim=c(0,10),xlim=c(min(sd$BP),max(sd$BP)),xlab='POS(hg19)',ylab='-log10(P)')
		abline(h=7.3,col='dark gray')
		  
		size<-input$bpEnd-input$bpStart
		par(mar=c(2,10,0,2),font.axis=1.5,font.lab=1.5)	
		plot(size,2,xlim=c(0,size),ylim=c(0,7),ann='F',axes='F',type='n',col='white')
		text(c(-3,-3,-3),c(2.75),labels=c('RefGene'),col=c('black'),cex=1.5,font=4,srt = 0,adj=1.85,xpd=T)	
		##refgene
		refgene=datasets[["ref"]]
		subref<-subset(refgene,(refgene[,1]==paste('chr',input$chr,sep='')&(refgene[,2]>=input$bpStart&refgene[,2]<=input$bpEnd))|(refgene[,1]==paste('chr',input$chr,sep='')&(refgene[,3]>=input$bpStart&refgene[,3]<=input$bpEnd)))
		d<-which(duplicated(subref[,6])=='TRUE')
		if(length(d)!=0){subref[d,6]<-''}
		if(nrow(subref)!=0){
		geneout<-NULL
		#gene figure
		for(i in 1:nrow(subref)){
		tmpstart<-as.character(subref[i,4])
		tmpend<-as.character(subref[i,5])
		start<-unlist(strsplit(tmpstart,',',fixed=T))
		end<-unlist(strsplit(tmpend,',',fixed=T))
		  for(j in 1:length(start)){
			  start1<-as.numeric(start[j])
			  end1<-as.numeric(end[j])
			  rect(start1-input$bpStart,2.5,end1-input$bpStart,3.0,col='black',border='black')
			  }	
		lines(c(subref[i,2]-input$bpStart,subref[i,3]-input$bpStart),c(2.75,2.75),col='black',lwd=1.6)
		text((subref[i,2]-input$bpStart+subref[i,3]-input$bpStart)/2,1,labels=subref[i,6],col='black',cex=1.0,font=4,srt=90,xpd=T)	  
		if(is.na(subref[i,6])!='TRUE'){geneout<-paste(geneout,',',subref[i,6],sep='')}
			}
		geneout<-sub(',','',geneout)
		}else{geneout<-'NA'}
	}}, res = 96)
}

shinyApp(ui, server)