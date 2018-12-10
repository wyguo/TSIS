
packages = c("shiny", "zoo","ggplot2","shinyFiles")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# require(tools)
library(shiny)
library(ggplot2)
library(zoo)
library(shinyFiles)
######################################################################################################
##TSIS
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = '*.R')) {
    #if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    #if(trace) cat("/n")
  }
}
sourceDir(path = 'R',encoding = 'UTF-8')


######################################################################################################

######################################################################################################

data.size.max=100

message('Starting Shiny App...')

##shiny app
shinyApp(options = list(launch.browser=T),
         ui = navbarPage("Time-series isoform switch",
                         # tags$head(tags$script(src="google-analytics.js")),
                         ##Page 1
                         tabPanel("Manual",
                                  tags$head(includeScript("https://raw.githubusercontent.com/wyguo/TSIS/master/inst/google-analytics.js")),
                                  # tags$head(tags$script(src="https://raw.githubusercontent.com/wyguo/TSIS/master/inst/google-analytics.js")),
                                  HTML('User manual: <a href="https://github.com/wyguo/TSIS" target="_blank">https://github.com/wyguo/TSIS</a>'),
                                  br(),
                                  HTML('Example data are in two csv files, data.exp.csv and mapping.csv: <a href="https://github.com/wyguo/TSIS/tree/master/vignettes/data" target="_blank">https://github.com/wyguo/TSIS/tree/master/vignettes/data</a>'),
                                  hr(),
                                  HTML('<b>If you use TSIS in your work, please cite:</b>'),
                                  br(),
                                  HTML('Wenbin Guo, Cristiane P. G. Calixto, John W.S. Brown, Runxuan Zhang, "TSIS: an R package to infer alternative splicing isoform switches for time-series data", Bioinformatics, https://doi.org/10.1093/bioinformatics/btx411, 2017.'),
                                  br(),
                                  br(),
                                  HTML('<b>The paper for TSIS application in RNA-seq data from study of Arabidopsis in response to cold:</b>'),
                                  br(),
                                  HTML('Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H., Zhang,R., and Brown,J.W.S. (2018) Rapid and dynamic alternative splicing impacts the Arabidopsis cold response transcriptome. Plant Cell.')
                                  
                         ),
                         ##Page 2
                         tabPanel("Isoform switch analysis",
                                  fluidPage(
                                    fluidRow(
                                      ##input gene expression data part
                                      column(3,
                                             titlePanel('Input data files'),
                                             wellPanel(
                                               
                                               ##input isoform expression data
                                               h4('Isoform expression data:'),
                                               # fileInput('filedata','Select the expression data file',
                                               #           accept = c(
                                               #             'text/csv',
                                               #             'text/comma-separated-values',
                                               #             'text/tab-separated-values',
                                               #             'text/plain',
                                               #             '.csv',
                                               #             '.tsv'
                                               #           )),
                                               shinyFilesButton(id = 'filedata', label = 'Select csv file', title = 'Isoform expression data', multiple = F),
                                               p('The expression input is a data table with columns of samples and rows of isoforms.'),
                                               br(),
                                               
                                               ##isoforms mapping data
                                               h4('Isoforms mapping data:'),
                                               # fileInput('filetarget','Select the mapping data file',
                                               #           accept = c(
                                               #             'text/csv',
                                               #             'text/comma-separated-values',
                                               #             'text/tab-separated-values',
                                               #             'text/plain',
                                               #             '.csv',
                                               #             '.tsv'
                                               #           )),
                                               shinyFilesButton(id = 'filetarget', label = 'Select csv file', title = 'Gene-transcript mapping data', multiple = F),
                                               p('The mapping input is a data table with first column of genes and second column of isoforms.')
                                             ),
                                             wellPanel(
                                               ##input subset of isoforms for investigation
                                               h4('Subset of isoforms for investigation:'),
                                               # fileInput('file.subtarget','Select the subset of isoforms data file',
                                               #           accept = c(
                                               #             'text/csv',
                                               #             'text/comma-separated-values',
                                               #             'text/tab-separated-values',
                                               #             'text/plain',
                                               #             '.csv',
                                               #             '.tsv'
                                               #           )),
                                               shinyFilesButton(id = 'filesubtarget', label = 'Select csv file', title = 'Select the subset of isoforms data file', multiple = F),
                                               HTML('Only the results of provided subset of isoforms and their isoform partners will be shown in the results.')
                                               
                                             ),
                                             titlePanel('Density/Frequency of switch'),
                                             wellPanel(
                                               plotOutput('density',height = "300px"),
                                               HTML('<b>Figure:</b> Density/Frequency plot of switch time points. The plot is made based
                                                    on the occurring time points of isoform switches after scoring and filtering.'),
                                               fluidRow(
                                                 column(6,
                                                        radioButtons("densityplot.type", label = h5("Select plot type:"),
                                                                     choices = list("Frequency" = 'frequency',"Density" = 'density'),
                                                                     selected = 'frequency',inline = F)
                                                 ),
                                                 column(6,
                                                        radioButtons("show.density.line", label = h5("Density/frequency line:"),
                                                                     choices = list("FALSE" = 'FALSE',"TRUE" = 'TRUE'),
                                                                     selected = 'FALSE',inline = F)
                                                 )
                                               ),
                                               fluidRow(
                                                 column(6,
                                                        radioButtons("densityplot.format", label = h5("Select format to save:"),
                                                                     choices = list("png" = "png", "pdf" = 'pdf'),
                                                                     selected = 'png',inline = T)
                                                 ),
                                                 column(6,
                                                        br(),
                                                        br(),
                                                        downloadButton('download.densityplot', 'Save',class="btn btn-primary")
                                                 )
                                               )
                                               )
                                    ),
                                    ##input taret part
                                    column(9,
                                           titlePanel('Parameter settings'),
                                           # tableOutput(outputId = 'test.table'),
                                           # tableOutput(outputId = 'test.table2'),
                                           wellPanel(
                                             fluidRow(
                                               h4('Scoring parameters'),
                                               # column(3,
                                               #        numericInput('t.start',label='Start time:',value=1)
                                               # ),
                                               # column(3,
                                               #        numericInput('t.end',label='End time:',value=26)
                                               # ),
                                               # column(3,
                                               #        numericInput('nrep',label='Replicates:',value=9)
                                               # ),
                                               column(3,
                                                      selectInput('method.intersection','Search intersections',c('Mean','Spline'))
                                               ),
                                               column(3,
                                                      numericInput('spline.df',label='Spline degree:',value=18)
                                               )
                                               
                                             ),
                                             
                                             actionButton('scoring','Scoring',icon("send outline icon"),class="btn btn-primary"),
                                             
                                             br(),
                                             br(),
                                             HTML('Press Scoring button to implement the scoring of isoform switches. The details of parameters:
                                                  <ul><li><b>Method for intersections:</b> Using either mean values or natural spline fitted smooth curves (see detail in ns() function in splines R package) of time-series expression
                                                  to determine the intersection points of isoforms.</li>
                                                  <li><b>Degree of spline:</b> The degree of spline in splines::ns() function.</li>
                                                  </ul>')
                                             ),
                                           wellPanel(
                                             fluidRow(
                                               h4('Filtering parameters'),
                                               column(3,
                                                      numericInput('prob.cutoff',label='Probability cutoff:',value=0.5)
                                               ),
                                               column(3,
                                                      numericInput('diff.cutoff',label='Difference cutoff:',value=1)
                                               ),
                                               column(3,
                                                      numericInput('pval.cutoff',label='p-value cutoff:',value=0.001)
                                               ),
                                               column(3,
                                                      numericInput('t.points.cutoff',label='Min time in interval:',value=2)
                                               ),
                                               column(3,
                                                      numericInput('cor.cutoff',label='Correlation cutoff:',value=0)
                                               ),
                                               column(3,
                                                      numericInput('x.lower.boundary',label='Lower time:',value=1)
                                               ),
                                               column(3,
                                                      numericInput('x.upper.boundary',label='Upper time:',value=26)
                                               ),
                                               column(3,
                                                      selectInput('sub.isoforms.ft','Subset of isoforms:',c('FALSE','TRUE'))
                                               ),
                                               column(3,
                                                      selectInput('max.ratio','Most abundant isoforms only',c('FALSE','TRUE'))
                                               )
                                               
                                             ),
                                             br(),
                                             actionButton('filtering','Filtering',icon("send outline icon"),class="btn btn-primary"),
                                             br(),
                                             br(),
                                             HTML('Press Filtering button to filter the scores. The details of parameters:
                                                  <ul><li><b>Probability cutoff:</b> The isoform switch probability/frequency cut-off for the column "prob" in the output table. </li>
                                                  <li><b>Difference cutoff:</b> The isoform switch difference cut-off for the column "diff" in the output table. </li>
                                                  <li><b>P-value cutoff:</b> The p-value cut-off of both columns "before.pval" and "after.pval" in the output table. </li>
                                                  <li><b>Min time in interval:</b> The minimum time points for both columns "before.t.points" and "after.t.points" in the output table.</li>
                                                  <li><b>Correlation cutoff:</b> The cut-off for Pearson correlation of isoform pairs.</li>
                                                  <li><b>Lower boundary of time, Upper boundary of time:</b> Specifies the time frame of interest to investigate the isoform switch. </li>
                                                  <li><b>Subset of isoforms:</b> Logical, output subset of the results only according to the subset isoform list in the input if selected. </li>
                                                  <li><b>Most abundant isoforms only:</b> Logical, only output the results of most abundant isoforms if selected. </li>
                                                  </ul>')
                                             ),
                                           titlePanel('Output results of isoform switch'),
                                           wellPanel(
                                             
                                             HTML('The following table shows the feature scores of isoform switches. The columns in the output table:
                                                  <ul><li><b>iso1, iso2:</b> the isoform pairs. </li>
                                                  <li><b>iso1.mean.ratio, iso2.mean.ratio:</b> The mean ratios of isoforms to their gene. </li>
                                                  <li><b>before.interval, after.interval:</b> The intervals before and after switch points. </li>
                                                  <li><b>x.value, y.value:</b> The values of x axis (time) and y axis (expression) coordinates of the switch points.</li>
                                                  <li><b>prob:</b> The probability/frequency of switch.</li>
                                                  <li><b>diff:</b> The sum of average sample differences before and after switch.</li>
                                                  <li><b>before.pval, after.pval:</b> The paired t-test p-values of the samples in the intervals before and after switch points.</li>
                                                  <li><b>before.t.points, after.t.points:</b> The number of time points in intervals before and after the switch points.</li>
                                                  <li><b>cor:</b> The Pearson correlation of iso1 and iso2.</li>
                                                  </ul>')
                                             )
                                             )
                                             ),
                                    
                                    fluidRow(
                                      column(12,
                                             wellPanel(
                                               HTML('<b>Table:</b> The TSIS analysis results'),
                                               br(),
                                               shiny::dataTableOutput('score.table')
                                             ))
                                    ),
                                    fluidRow(
                                      column(10,
                                             div(align='right',downloadButton('download.scores', 'Download score table',class="btn btn-primary"))
                                      ),
                                      column(2,
                                             div(align='right',downloadButton('download.genes', 'Download gene names',class="btn btn-primary"))
                                      )
                                    )
                                             )
                                  
                                           ),
                         ##Page 3
                         tabPanel("Visualization",
                                  fluidRow(
                                    titlePanel('Switch plots'),
                                    column(2,
                                           wellPanel(
                                             p(h4('Input isoform names:')),
                                             textInput('iso1', label = 'Isoform 1: ', value = ''),
                                             textInput('iso2', label = 'Isoform 2: ', value = ''),
                                             numericInput('prob.cutoff.switch.points',label = 'Show switch poitns with prob >:',value = 0.5),
                                             selectInput('ribbon.plot','Plot types',c('Error bar','Ribbon')),
                                             radioButtons("show.scores", label = h4("Show feature labels:"),
                                                          choices = list("TRUE" = 'TRUE', "FALSE" = "FALSE"),
                                                          selected = 'TRUE',inline = T),
                                             radioButtons("show.color.region", label = h4("Show region:"),
                                                          choices = list("TRUE" = 'TRUE', "FALSE" = "FALSE"),
                                                          selected = 'TRUE',inline = T),
                                             actionButton('plot1iso','Plot',icon("send outline icon"),class="btn btn-primary"),
                                             br(),
                                             br(),
                                             # p(h4('Select format to save:')),
                                             # selectInput('plot1.format','',c('html','png','pdf')),
                                             radioButtons("plot1.format", label = h4("Select format to save:"),
                                                          choices = list("png" = "png", "pdf" = 'pdf'),
                                                          selected = 'png',inline = T),
                                             downloadButton('download.1plot', 'Save',class="btn btn-primary")
                                             
                                           )
                                    ),
                                    column(10,
                                           wellPanel(
                                             plotOutput('plot.isoform',width = 'auto',height = '550px')
                                           )
                                    )
                                  ),
                                  fluidRow(
                                    titlePanel('Isoform switch plots in batch'),
                                    
                                    column(2,
                                           wellPanel(
                                             p(h4('Isoforms pairs to plot:')),
                                             numericInput('topn','Top n isofomrs:',value = 10),
                                             # selectInput('plotmulti.format','Select format to save:',c('png','pdf'))
                                             radioButtons("plotmulti.format", label = h4("Select format to save:"),
                                                          choices = list("png" = "png", "pdf" = 'pdf'),
                                                          selected = 'png',inline = T)
                                             
                                           )),
                                    column(5,
                                           wellPanel(
                                             textInput('folder2save', label = 'Folder to save: ', value = 'figure',width='400px'),
                                             HTML('Note: typing a folder name to save the multiple plots. If the folder does not exist in the working directory, a
                                                  new folder is created with the name.'),
                                             br(),
                                             br(),
                                             actionButton('plotmultiiso','Plot',icon("send outline icon"),buttonType ='button',class="btn btn-primary"),
                                             h5(textOutput("directorypath", container = span))
                                             
                                             
                                             ))),
                                  fluidRow(
                                    column(12,
                                           wellPanel(
                                             h4('Isoform pairs to plot'),
                                             shiny::dataTableOutput('isoform.pairs.to.plot'),
                                             textOutput('isoform.pairs.to.plot.loops')
                                           )
                                    )
                                  )
                         )
                                           ),
         
         
         
         ##server function
         
         
         
         server = function(input, output,session) {
           options(shiny.maxRequestSize = data.size.max*1024^2)
           shinyFileChoose(input, 'filedata', roots=c(wd='/srv/shiny-server'), filetypes=c('txt', 'csv'))
           shinyFileChoose(input, 'filetarget', roots=c(wd='/srv/shiny-server'), filetypes=c('txt', 'csv'))
           shinyFileChoose(input, 'filesubtarget', roots=c(wd='/srv/shiny-server'), filetypes=c('txt', 'csv'))
           
           data.exp0<-reactive({
             inFile.data <- parseFilePaths(roots=c(wd='/srv/shiny-server'), input$filedata)
             if (is.null(inFile.data))
               return(NULL)
             data.exp0<-read.csv(file=as.character(inFile.data$datapath),header=F)
             return(data.exp0)
             # save(data.exp0,file='data/data.exp0.RData')
           })
           
    
           times<-reactive({
             if (is.null(data.exp0))
               return(NULL)
             
             times<-as.numeric(as.vector(t(data.exp0()[2,-1])))
             return(times)
             
           })
           
           
           data.exp<-reactive({
             if (is.null(data.exp0()))
               return(NULL)
             
             data.exp<-data.exp0()
             values<-data.exp[-c(1:2),-1]
             values<-data.frame(apply(values,2,as.numeric))
             rownames(values)<-as.character(data.exp[-c(1:2),1])
             colnames(values)<-paste0(as.character(t(data.exp[1,-1])),'_',as.character((t(data.exp[2,-1]))))
             data.exp<-values
             data.exp<-na.omit(data.exp)
             return(data.exp)
             
           })
           
           mapping<-reactive({
             inFile.mapping <- parseFilePaths(roots=c(wd='/srv/shiny-server'), input$filetarget)
             if (is.null(inFile.mapping))
               return(NULL)
             mapping <- read.csv(file=as.character(inFile.mapping$datapath),header=T)
             mapping<-na.omit(mapping)
             return(mapping)
           })
           # 
           # output$test.table <- renderTable({
           #   data.frame(mapping()[1:10,])
           # })
           # 
           # output$test.table2 <- renderTable({
           #   data.frame(data.exp()[1:10,1:3])
           # })
           # 
           
           sub.isoform.list<-reactive({
             # infile.subtarget<-input$file.subtarget
             inFile.sublist <- parseFilePaths(roots=c(wd='/srv/shiny-server'), input$filesubtarget)
             if (is.null(inFile.sublist))
               return(NULL)
             if(input$sub.isoforms.ft=='TRUE' & !is.null(inFile.sublist))
               sub.isoform.list<-as.vector(t(read.csv(file=as.character(inFile.sublist$datapath),header=T))) else sub.isoform.list<-NULL
               
           })
           
           
           
           
           scores<-eventReactive(input$scoring,{
             if (is.null(data.exp()) | is.null(mapping()))
               return(NULL)
             
             ##parameter for itch()
             # t.start<-input$t.start
             # t.end<-input$t.end
             # nrep<-input$nrep
             min.t.points<-input$t.points.cutoff
             min.difference<-input$diff.cutoff
             ##parameters for
             scores<-iso.switch.shiny(data.exp=data.exp(),mapping=mapping(),times = times(),
                                      min.t.points = min.t.points,min.difference = min.difference,rank = F,
                                      spline = input$method.intersection=='Spline',spline.df = input$spline.df)
             
           })
           
           scores.filtered<-eventReactive(input$filtering,{
             if (is.null(data.exp()) | is.null(mapping()))
               return()
             
             t.points.cutoff<-input$t.points.cutoff
             prob.cutoff<-input$prob.cutoff
             diff.cutoff<-input$diff.cutoff
             pval.cutoff<-input$pval.cutoff
             cor.cutoff<-input$cor.cutoff
             x.value.limit<-c(input$x.lower.boundary,input$x.upper.boundary)
             #
             scores.filtered<-score.filter(scores = scores(),t.points.cutoff=t.points.cutoff,prob.cutoff=prob.cutoff,diff.cutoff=diff.cutoff,
                                           pval.cutoff=pval.cutoff,cor.cutoff = cor.cutoff,x.value.limit=x.value.limit,sub.isoform.list=sub.isoform.list(),
                                           sub.isoform=(input$sub.isoforms.ft=='TRUE'),max.ratio=(input$max.ratio=='TRUE'),
                                           data.exp=data.exp(),mapping=mapping()
             )
             
             
             
             
             
           })
           
           
           score.show<-reactiveValues(scores=NULL)
           
           observeEvent(input$scoring, {
             if (is.null(data.exp()) | is.null(mapping()))
               return()
             
             x<- scores()
             if(nrow(x)>1)
               x[,-c(1,2,5,6)]<-apply(x[,-c(1,2,5,6)],2,function(x) as.numeric(format(x,digits = 3)))
             
             rownames(x)<-NULL
             score.show$scores<-x[,c('iso1','iso2','iso1.mean.ratio','iso2.mean.ratio',
                                     'before.interval','after.interval','x.value','y.value','prob','diff','before.pval','after.pval','before.t.points','after.t.points','cor')]
             
             
           })
           
           
           observeEvent(input$filtering, {
             if (is.null(data.exp()) | is.null(mapping()))
               return()
             
             x<- scores.filtered()
             if(nrow(x)>1)
               x[,-c(1,2,5,6)]<-apply(x[,-c(1,2,5,6)],2,function(x) as.numeric(format(x,digits = 3)))
             
             rownames(x)<-NULL
             score.show$scores<-x[,c('iso1','iso2','iso1.mean.ratio','iso2.mean.ratio',
                                     'before.interval','after.interval','x.value','y.value','prob','diff','before.pval','after.pval','before.t.points','after.t.points','cor')]
             
             
           })
           
           
           
           
           output$score.table <- shiny::renderDataTable(options=list(pageLength=20,
                                                                     aoColumnDefs = list(list(sClass="alignLeft",aTargets="_all")),
                                                                     scrollX=TRUE,scrollY="500px"),{
                                                                       if(is.null(score.show$scores))
                                                                         return()
                                                                       
                                                                       score.show$scores
                                                                       
                                                                     })
           
           
           
           ##save scores to csv files
           ##save p-value table
           output$download.scores <- downloadHandler(
             filename=function(){
               'scores.csv'
             },
             content=function(file){
               write.csv(score.show$scores,file,row.names = F)
             }
           )
           
           output$download.genes <- downloadHandler(
             filename=function(){
               'genes.csv'
             },
             content=function(file){
               write.csv(data.frame(genes=unique(mapping()[which(as.vector(mapping()[,2]) %in% unique(c(as.vector(score.show$scores$iso1),as.vector(score.show$scores$iso2)))),1])),
                         file,row.names = F)
             }
           )
           
           ##plot the density
           # height = 400, width = 600
           output$density <- renderPlot({
             if(is.null(score.show$scores))
               return()
             
             switch.density(x=score.show$scores$x.value,time.points = unique(times()),
                            plot.type = input$densityplot.type,make.plotly = F,
                            show.line = input$show.density.line,title = '',autosize = F,width=250,height=250)
           })
           
           
           output$download.densityplot <- downloadHandler(
             
             filename = function() {
               paste0('Density plot.',input$densityplot.format)
             },
             content = function(file,format=input$densityplot.format) {
               ggsave(file,
                           switch.density(x=score.show$scores$x.value,time.points = unique(times()),
                                          plot.type = input$densityplot.type,make.plotly = F,
                                          show.line = input$show.density.line,title = ''),
                           width = 16,height = 12,units = "cm")
             })
           
           
           #########
           
           
           g<-eventReactive(input$plot1iso,{
             if(is.null(input$iso1) | is.null(input$iso2) | is.null(data.exp()))
               return()
             
             iso1<-trimws(input$iso1)
             iso2<-trimws(input$iso2)
             g<-plotTSIS(data2plot = data.exp()[c(iso1,iso2),],
                         iso1 = iso1,make.plotly = F,
                         show.region = input$show.color.region,
                         iso2 = iso2,
                         scores = scores.filtered(),
                         show.scores = input$show.scores,times = times(),
                         x.lower.boundary=input$x.lower.boundary,
                         x.upper.boundary=input$x.upper.boundary,
                         prob.cutoff=input$prob.cutoff.switch.points,
                         y.lab = 'Expression',spline=(input$method.intersection=='Spline'),spline.df = input$spline.df,ribbon.plot = (input$ribbon.plot=='Ribbon')
                         
             )
           })
           
           output$plot.isoform<-renderPlot({
             if(is.null(g()))
               return()
             
             ggplotly(g(),autosize = F,width=1000,height=500)
           })
           #
           
           output$download.1plot <- downloadHandler(
             
             filename = function() {
               paste0(input$iso1,' vs ',input$iso2,'.',input$plot1.format)
             },
             content = function(file,format=input$plot1.format) {
               ggsave(file, g(),width = 25,height = 12,units = "cm")
             })
           
           
           topn<-eventReactive(input$plotmultiiso,{
             topn<-input$topn
           })
           
           plotmulti.format<-eventReactive(input$plotmultiiso,{
             plotmulti.format<-input$plotmulti.format
           })
           
           
           
           folderInput <- eventReactive(input$plotmultiiso,{
             x<-paste0('/srv/shiny-server','/',trimws(input$folder2save))
             if(!file.exists(x))
               dir.create(x)
             x
           })
           
           
           output$directorypath = renderText({
             paste0('Plots are saved in: "', folderInput(),'"')
           })
           
           
           
           data2plot<-eventReactive(input$plotmultiiso,{
             if(is.null(scores.filtered()))
               return()
             
             topn<-topn()
             topn<-min(topn,nrow(scores.filtered()))
             
             data2plot<-score.show$scores[1:topn,]
             
           })
           
           
           
           
           output$isoform.pairs.to.plot <- shiny::renderDataTable(options=list(pageLength=10,
                                                                               aoColumnDefs = list(list(sClass="alignLeft",aTargets="_all")),
                                                                               scrollX=TRUE,scrollY="400px"),{
                                                                                 if(is.null(score.show$scores))
                                                                                   return()
                                                                                 
                                                                                 score.show$scores[1:input$topn,]
                                                                                 
                                                                               })
           
           
           output$isoform.pairs.to.plot.loops<-renderText({
             if(is.null(data2plot()))
               return()
             
             start.time <- Sys.time()
             x<-data2plot()
             
             
             iso1s<-as.character(x[,1])
             iso2s<-as.character(x[,2])
             
             # c(iso1s,iso2s)
             withProgress(message = 'Making plots...',value=0,{
               for(i in 1:length(iso1s)){
                 gs<-plotTSIS(data2plot = data.exp()[c(iso1s[i],iso2s[i]),],line.width= 1,
                              iso1 = NULL,make.plotly = F,
                              iso2 = NULL,
                              scores = data2plot(),times=times(),
                              x.lower.boundary=input$x.lower.boundary,
                              x.upper.boundary=input$x.upper.boundary,
                              prob.cutoff=input$prob.cutoff.switch.points,
                              y.lab = 'Expression',spline=(input$method.intersection=='Spline'),spline.df = input$spline.df,ribbon.plot = (input$ribbon.plot=='Ribbon')
                 )
                 
                 
                 plot.name<-paste0(folderInput(),'/',iso1s[i],' vs ',iso2s[i],'.',plotmulti.format())
                 ggsave(plot.name,gs,width = 25,height = 12,units = "cm")
                 
                 # unlink(plot.name)
                 
                 incProgress(1/length(iso1s), detail = paste(i, ' of ', length(iso1s)))
                 Sys.sleep(0)
                 
               }
               
             })
             
             end.time <- Sys.time()
             time.taken <- end.time - start.time
             paste0('Done. Time for ploting: ',format(time.taken,digits=4,unit='auto'))
             
           })
         })


