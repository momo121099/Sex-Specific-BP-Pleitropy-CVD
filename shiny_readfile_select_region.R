library(shiny)
library(data.table)  # for fread, which is generally faster and can handle larger files

readData <- function() {
  data.sbp <- fread("SBP_female_male_bolt_imputed_result_MAF0.01_short.txt", header = TRUE)
  data.dbp <- fread("DBP_female_male_bolt_imputed_result_MAF0.01_short.txt", header = TRUE)
  data.pp <- fread("PP_female_male_bolt_imputed_result_MAF0.01_short.txt", header = TRUE)
  return(list(data.sbp = data.sbp, data.dbp = data.dbp, data.pp = data.pp))
}

# Create a global variable to store the data
data <- readData()

ui <- fluidPage(
  titlePanel("Subset Data by CHR and BP Range"),
  sidebarLayout(
    sidebarPanel(
      selectInput("select", h3("BP trait"), choices = list("SBP" = 'SBP', "DBP" = 'DBP',"PP" = 'PP'), selected = 'SBP'),
      numericInput("chr", "Select CHR:", 13),
      numericInput("bpStart", "Select BP Start:", 110546007),
      numericInput("bpEnd", "Select BP End:", value = 111046007)
    ),
    mainPanel(
      plotOutput("plot", click = "plot_click", height = "800px"),
      div(style = "width: 100%; text-align: center;",
          h3("Top SNP in male", style = "font-size: 18px;"),
          div(style = "display: inline-block; margin-left: auto; margin-right: auto;",
              tableOutput("subsetTable1")
          )  # Closing div for the male table
      ),  # Closing div for the male section
      div(style = "width: 100%; text-align: center;",
          h4("Top SNP in female", style = "font-size: 18px;"),
          div(style = "display: inline-block; margin-left: auto; margin-right: auto;", 
              tableOutput("subsetTable2")
          )  # Closing div for the female table
      )  # Closing div for the female section
    )
  )
)


server <- function(input, output) {
  # Reactive expression for the subset data
  subset_data <- reactive({
    # Ensure that the input values are numeric and within the expected range
    chr <- as.numeric(input$chr)
    bpStart <- as.numeric(input$bpStart)
    bpEnd <- as.numeric(input$bpEnd)
    
    # Filter the data based on the input values
    selected_data <- switch(input$select,
      SBP = data$data.sbp,
      DBP = data$data.dbp,
      PP = data$data.pp
    )
    
    selected_data[CHR == chr & BP >= bpStart & BP <= bpEnd, ]
  })
  
  # Output the subset data as a table
  output$subsetTable1 <- renderTable({
    tb=subset_data()
    tb1=tb[order(tb$P_BOLT_LMM.male)[1],]
    tb1$P_BOLT_LMM.male=format(tb1$P_BOLT_LMM.male, scientific = TRUE)
    tb11=data.frame(tb1$CHR,tb1$BP,tb1$A1FREQ.male,tb1$BETA.male,tb1$SE.male,tb1$P_BOLT_LMM.male)
    colnames(tb11)=c("CHR","BP","FREQ","BETA","SE","P")
    tb11

  })

  output$subsetTable2 <- renderTable({
    tb=subset_data()
    tb2=tb[order(tb$P_BOLT_LMM.female)[1],]
    tb2$P_BOLT_LMM.female=format(tb2$P_BOLT_LMM.female, scientific = TRUE)
    tb22=data.frame(tb2$CHR,tb2$BP,tb2$A1FREQ.female,tb2$BETA.female,tb2$SE.female,tb2$P_BOLT_LMM.female)
    colnames(tb22)=c("CHR","BP","FREQ","BETA","SE","P")
    tb22
  })

  # Output the plot
  output$plot <- renderPlot({
    # Access the reactive subset_data within the renderPlot
    sd <- subset_data()
    
    # Check if the columns exist before plotting
    if("P_BOLT_LMM.male" %in% names(sd) && "P_BOLT_LMM.female" %in% names(sd)) {
      par(mfrow=c(2,1))
      max=max(c(-log10(sd$P_BOLT_LMM.male),-log10(sd$P_BOLT_LMM.female)))
      plot(sd$BP, -log10(sd$P_BOLT_LMM.male), main = "Male-only Genetic Association",pch=19,col='blue',ylim=c(0,max),xlim=c(min(sd$BP),max(sd$BP)),xlab='POS(hg19)',ylab='-log10(P)')
      plot(sd$BP, -log10(sd$P_BOLT_LMM.female), main = "Female-only Genetic Association",pch=19,col='purple',ylim=c(0,max),xlim=c(min(sd$BP),max(sd$BP)),xlab='POS(hg19)',ylab='-log10(P)')
    }
  }, res = 96)
}

shinyApp(ui, server)
