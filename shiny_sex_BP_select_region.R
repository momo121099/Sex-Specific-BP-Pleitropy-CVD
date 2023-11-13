library(shiny)
library(data.table)  # for fread, which is generally faster and can handle larger files

readData <- function() {

  sbp.f<-fread('UKB_SBP_female_bolt_imputed_result_MAF0.01.tsv',header=T)
  sbp.m<-fread('UKB_SBP_male_bolt_imputed_result_MAF0.01.tsv',header=T)
  dbp.f<-fread('UKB_DBP_female_bolt_imputed_result_MAF0.01.tsv',header=T)
  dbp.m<-fread('UKB_DBP_male_bolt_imputed_result_MAF0.01.tsv',header=T)
  pp.f<-fread('UKB_PP_female_bolt_imputed_result_MAF0.01.tsv',header=T)
  pp.m<-fread('UKB_PP_male_bolt_imputed_result_MAF0.01.tsv',header=T)

  data.sbp <- merge(sbp.f, sbp.m, by="pos.hg19", all=FALSE, suffixes=c("_female","_male"))
  data.dbp <- merge(dbp.f, dbp.m, by="pos.hg19", all=FALSE, suffixes=c("_female","_male"))
  data.pp <- merge(pp.f, pp.m, by="pos.hg19", all=FALSE, suffixes=c("_female","_male")) 
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
    
    selected_data[chromosome_female == chr & base_pair_location_female >= bpStart & base_pair_location_female <= bpEnd, ]
  })
  
  # Output the subset data as a table
  output$subsetTable1 <- renderTable({
    tb=subset_data()
    tb1=tb[order(tb$p_value_male)[1],]
    tb1$p_value_male=format(tb1$p_value_male, scientific = TRUE)
    tb11=data.frame(tb1$chromosome_male,tb1$base_pair_location_male,tb1$effect_allele_male,tb1$beta_male,tb1$standard_error_male,tb1$p_value_male)
    colnames(tb11)=c("CHR","BP","FREQ","BETA","SE","P")
    tb11

  })

  output$subsetTable2 <- renderTable({
    tb=subset_data()
    tb2=tb[order(tb$p_value_female)[1],]
    tb2$p_value_female=format(tb2$p_value_female, scientific = TRUE)
    tb22=data.frame(tb2$chromosome_female,tb2$base_pair_location_female,tb2$effect_allele_female,tb2$beta_female,tb2$standard_error_female,tb2$p_value_female)
    colnames(tb22)=c("CHR","BP","FREQ","BETA","SE","P")
    tb22
  })

  # Output the plot
  output$plot <- renderPlot({
    # Access the reactive subset_data within the renderPlot
    sd <- subset_data()
    
    # Check if the columns exist before plotting
    if("p_value_male" %in% names(sd) && "p_value_female" %in% names(sd)) {
      par(mfrow=c(2,1))
      max=max(c(-log10(sd$p_value_male),-log10(sd$p_value_female)))
      plot(sd$base_pair_location_female, -log10(sd$p_value_male), main = "Male-only Genetic Association",pch=19,col='blue',ylim=c(0,max),xlim=c(min(sd$base_pair_location_female),max(sd$base_pair_location_female)),xlab='POS(hg19)',ylab='-log10(P)')
      plot(sd$base_pair_location_female, -log10(sd$p_value_female), main = "Female-only Genetic Association",pch=19,col='purple',ylim=c(0,max),xlim=c(min(sd$base_pair_location_female),max(sd$base_pair_location_female)),xlab='POS(hg19)',ylab='-log10(P)')
    }
  }, res = 96)
}

shinyApp(ui, server)