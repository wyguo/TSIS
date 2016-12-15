#' Isoform switch analysis and visualization with Shiny App
#'
#' @param data.size.max maximum size limit for unload files in Shiny. Default is 100 (MB).
#'
#' @return the Shiny App.
#'
#' @examples TSIS.app()
#'
#' @seealso \code{\link{shiny}}
#'
#' @export


TSIS.app <- function(data.size.max=100) {
  # require(tools)
  require(shiny)
  # require(shinyFiles)
  library(shinythemes)
  library(plotly)
  # library(TSIS)



  shinyApp(options = list(launch.browser=T),
           ui = navbarPage("Time-series isoform switch",

                           ##Page 1
                           tabPanel("Manual",
                                    htmlOutput("tutorial")

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
                                                         fileInput('filedata','Select the expression data file',
                                                                   accept = c(
                                                                     'text/csv',
                                                                     'text/comma-separated-values',
                                                                     'text/tab-separated-values',
                                                                     'text/plain',
                                                                     '.csv',
                                                                     '.tsv'
                                                                   )),
                                                         p('The expression input is a data frame with columns of samples and rows are isoforms.'),
                                                         br(),

                                                         ##isoforms mapping data
                                                         h4('Isoforms mapping data:'),
                                                         fileInput('filetarget','Select the mapping data file',
                                                                   accept = c(
                                                                     'text/csv',
                                                                     'text/comma-separated-values',
                                                                     'text/tab-separated-values',
                                                                     'text/plain',
                                                                     '.csv',
                                                                     '.tsv'
                                                                   )),
                                                         p('The mapping input is a data frame with first column of genes and second column of its isoforms.')
                                                       ),
                                                       wellPanel(
                                                         ##input subset of isoforms for investigation
                                                         h4('Subset of isoforms for investigation:'),
                                                         fileInput('file.subtarget','Subset of isoforms',
                                                                   accept = c(
                                                                     'text/csv',
                                                                     'text/comma-separated-values',
                                                                     'text/tab-separated-values',
                                                                     'text/plain',
                                                                     '.csv',
                                                                     '.tsv'
                                                                   )),
                                                         HTML('Only the results of provided subset of isoforms and their isoform partners will be shown in the results.'),
                                                         br(),
                                                         HTML('<b>Note:</b> The input data must be in *.csv format for loading convenience. The details of how to get example datasets are in
                                                               the manual section.' )

                                                       ),
                                                       wellPanel(
                                                         plotlyOutput('density',height = "300px"),
                                                         HTML('<b>Figure:</b> Density plot of switch time points. The plot is made based
                                                              on the x.value in the scores output table, i.e. the occurring time points of isoform switches.')
                                                       )
                                                ),
                                                ##input taret part
                                                column(9,
                                                       titlePanel('Parameter settings'),
                                                       wellPanel(
                                                         fluidRow(
                                                           h4('Parameters to calculate TS isoform switch scores'),
                                                           column(3,
                                                                  numericInput('t.start',label='Start time point:',value=1)
                                                           ),
                                                           column(3,
                                                                  numericInput('t.end',label='End time point:',value=26)
                                                           ),
                                                           column(3,
                                                                  numericInput('nrep',label='Number of replicates:',value=9)
                                                           ),
                                                           column(3,
                                                                  numericInput('min.t.points',label='Min time points of interval:',value=2)
                                                           ),
                                                           column(3,
                                                                  numericInput('min.distance',label='Min distance of isofomrs:',value=1)
                                                           ),
                                                           column(3,
                                                                  selectInput('rank','Using rank of expression',c('FALSE','TRUE'))
                                                           ),
                                                           column(3,
                                                                  selectInput('method.intersection','Method for intersections',c('Mean','Spline'))
                                                           ),
                                                           column(3,
                                                                  numericInput('spline.df',label='Degree of spline:',value=18)
                                                           )

                                                         ),

                                                         actionButton('scoring','Scoring',icon("send outline icon"),class="btn btn-primary"),

                                                         br(),
                                                         br(),
                                                         HTML('Press Scoring buttom to implement the scoring of isoform switches. The details of parameters:
                                                              <ul><li><b>Start time point and end time point:</b> The start time point and end time point of the time course. The time steps are assumed to be 1. </li>
                                                              <li><b>Number of replicates:</b> The number of replicates for the time points. </li>
                                                              <li><b>Min time points of interval:</b> Pre-filter the isoform pairs if all the interval only with time points <
                                                              Min time points of interval. </li>
                                                              <li><b>Min distance of isoforms:</b> Pre-filter the isoform pairs if the average distance between two isoforms at
                                                              every sample point < Min distance of isoforms.</li>
                                                              <li><b>Using rank of expression:</b> Logical, to take ranks of expression in each sample or not. </li>
                                                              <li><b>Method for intersections:</b> Using either mean values or natural spline fitted smooth curves (see detail of splines::ns() function in splines R package) of time-series expresssion
                                                              to determine the intersection poitns of time coures for the isoforms.</li>
                                                              <li><b>Degree of spline:</b> The degree of spline in splines::ns() function.</li>
                                                              </ul>'),
                                                         HTML('Note: the parameters for pre-filtering will reduce the computational intensity for large scale datasets.')
                                                         ),
                                                       wellPanel(
                                                         fluidRow(
                                                           h4('Parameters to filter TS isoform switch scores'),
                                                           column(3,
                                                                  numericInput('prob.cutoff',label='Probability cutoff:',value=0.5)
                                                           ),
                                                           column(3,
                                                                  numericInput('dist.cutoff',label='Distance cutoff:',value=1)
                                                           ),
                                                           column(3,
                                                                  numericInput('pval.cutoff',label='P-value cutoff:',value=0.01)
                                                           ),
                                                           column(3,
                                                                  numericInput('t.points.cutoff',label='Min time points of interval:',value=2)
                                                           ),
                                                           column(3,
                                                                  numericInput('cor.cutoff',label='Correlation cutoff:',value=0.5)
                                                           ),
                                                           column(3,
                                                                  numericInput('x.lower.boundary',label='Lower boundary of time:',value=9)
                                                           ),
                                                           column(3,
                                                                  numericInput('x.upper.boundary',label='Upper boundary of time:',value=17)
                                                           ),
                                                           column(3,
                                                                  selectInput('sub.isoforms.ft','Subset of isoforms:',c('FALSE','TRUE'))
                                                           ),
                                                           column(3,
                                                                  selectInput('max.ratio','Pairs of maximun isoform ratio:',c('FALSE','TRUE'))
                                                           )

                                                         ),
                                                         br(),
                                                         actionButton('filtering','Filtering',icon("send outline icon"),class="btn btn-primary"),
                                                         br(),
                                                         br(),
                                                         HTML('Press Filtering buttom to filter the scores. The details of parameters:
                                                              <ul><li><b>Probability cutoff:</b> The switch isoform probability/frequency cut-off for the column "prob" in the output table. </li>
                                                              <li><b>Distance cutoff:</b> The switch isoform distance cut-off for the column "dist" in the output table. </li>
                                                              <li><b>P-value cutoff:</b> The p-value cut-off for both columns "left.pval" and "right.pval" in the output table. </li>
                                                              <li><b>Min time points of interval:</b> The minimun time points for both columns "left.t.points" and "right.t.points" in the output table.</li>
                                                              <li><b>Correlation cutoff:</b> The cut-off for Pearson correaltion of isoform pairs.</li>
                                                              <li><b>Lower boundary of time, Upper boundary of time:</b> The lower and upper boundary of time duration for investigation. </li>
                                                              <li><b>Subset of isoforms:</b> Logical, to output only subset of the results or not? If TRUE, the subset of isoform list must be provided. </li>
                                                              <li><b>Pairs of maximun isoform ratio:</b> Logical, to output subset of the results of isoforms with maximum ratio or not? </li>
                                                              </ul>')
                                                         )
                                                       )
                                                       ),
                                              fluidRow(
                                                titlePanel('Input datasets'),
                                                column(9,
                                                       wellPanel(
                                                         h4('Partial expression data table'),
                                                         shiny::dataTableOutput('expression.table')
                                                       )
                                                ),
                                                column(3,
                                                       wellPanel(
                                                         h4('Mapping data table'),
                                                         shiny::dataTableOutput('target.table')
                                                       )
                                                )
                                              ),
                                              fluidRow(
                                                titlePanel('Output scores for isoform switch'),
                                                column(12,
                                                       wellPanel(
                                                         HTML('The columns in the output table:
                                                              <ul><li><b>iso1, iso2:</b> the isoform pairs. </li>
                                                              <li><b>iso1.mean.ratio, iso2.mean.ratio:</b> The ratios of isoforms to their gene. </li>
                                                              <li><b>left.interval, right.interval:</b>The breaks of the left and right intervals before and after swtich, respectively. </li>
                                                              <li><b>x.value, y.value:</b> The x axis and y axis value of the swtich points in the plot coordinates.</li>
                                                              <li><b>prob:</b> The probability of switch, which is defined as  <i>|probability(iso1>iso2|in left.interval)+probability(iso2>iso1|in right.interval)-1|</i>.</li>
                                                              <li><b>dist:</b> The average distance of switch, which is defined as  <i>|mean distance(iso1,iso2|in left.interval)|+|mean distance(iso1,iso2|in right.interval)|</i>.</li>
                                                              <li><b>left.pval, right.pval:</b> The paired t-test p-values for the samples in the left and rigth intervals for the switch points.</li>
                                                              <li><b>left.t.points, right.t.points:</b> The number of time points in the left and rigth intervals for the switch points.</li>
                                                              <li><b>cor:</b> The Pearson correlation of isoforms iso1 and iso2.</li>
                                                              </ul>'),

                                                         shiny::dataTableOutput('score.table')
                                                       ))
                                              ),
                                              fluidRow(
                                                column(12,
                                                       div(align='right',downloadButton('download.scores', 'Download'))
                                                )
                                              )
                                                       )

                                                       ),
                           ##Page 3
                           tabPanel("Switch visulisation",
                                    fluidRow(
                                      titlePanel('Switch plots'),
                                      column(2,
                                             wellPanel(
                                               p(h4('Input isoform names:')),
                                               textInput('iso1', label = 'Isoform 1: ', value = ''),
                                               textInput('iso2', label = 'Isoform 2: ', value = ''),
                                               numericInput('prob.cutoff.switch.points',label = 'Show switch poitns with prob >:',value = 0.5),
                                               selectInput('ribbon.plot','Plot types',c('Error bar','Ribbon')),
                                               radioButtons("show.scores", label = h4("Show score labels:"),
                                                            choices = list("TRUE" = 'TRUE', "FALSE" = "FALSE"),
                                                            selected = 'TRUE',inline = T),
                                               actionButton('plot1iso','Plot',icon("send outline icon"),class="btn btn-primary"),
                                               br(),
                                               br(),
                                               # p(h4('Select format to save:')),
                                               # selectInput('plot1.format','',c('html','png','pdf')),
                                               radioButtons("plot1.format", label = h4("Select format to save:"),
                                                            choices = list("html" = 'html', "png" = "png", "pdf" = 'pdf'),
                                                            selected = 'html',inline = T),
                                               downloadButton('download.1plot', 'Save',class="btn btn-primary")

                                             )
                                      ),
                                      column(10,
                                             wellPanel(
                                               plotlyOutput('plot.isoform',width = 'auto',height = '550px')
                                             )
                                      )
                                    ),
                                    fluidRow(
                                      titlePanel('Multiple plots for switch'),

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

             output$tutorial<-renderUI({
               includeHTML("vignettes/tutorial-shiny.html")
             })


             options(shiny.maxRequestSize = data.size.max*1024^2)



             data.exp<-reactive({
               infile.data<-input$filedata
               if (is.null(infile.data))
                 return(NULL)

               data.exp<-read.csv(file=infile.data$datapath,header=T)
               rownames(data.exp)<-data.exp[,1]
               data.exp<-data.exp[,-1]
             })

             mapping<-reactive({
               infile.target<-input$filetarget
               if (is.null(infile.target))
                 return(NULL)

               mapping<-read.csv(file=infile.target$datapath,header=T)
             })

             sub.isoform.list<-reactive({
               infile.subtarget<-input$file.subtarget

               if(input$sub.isoforms.ft=='TRUE' & !is.null(infile.subtarget))
                 sub.isoform.list<-as.vector(t(read.csv(file=infile.subtarget$datapath,header=T))) else sub.isoform.list<-NULL

             })


             output$expression.table <- shiny::renderDataTable(options=list(pageLength=10,aoColumnDefs = list(list(sClass="alignLeft",aTargets="_all")) ),{
               if (is.null(data.exp()))
                 return()

               if(ncol(data.exp())>5)
                 data.frame(isoforms=rownames(data.exp()),data.exp()[,1:4],`...`='...') else data.frame(isoforms=rownames(data.exp()),data.exp())
             })

             output$target.table <- shiny::renderDataTable(options=list(pageLength=10,aoColumnDefs = list(list(sClass="alignLeft",aTargets="_all")) ),{
               if (is.null(mapping()))
                 return()

               mapping()
             })




             scores<-eventReactive(input$scoring,{
               if (is.null(data.exp()) | is.null(mapping()))
                 return()

               ##parameter for itch()
               t.start<-input$t.start
               t.end<-input$t.end
               nrep<-input$nrep
               min.t.points<-input$min.t.points
               min.distance<-input$min.distance
               ##parameters for
               scores<-iso.switch.shiny(data.exp=data.exp(),mapping=mapping(),
                                        t.start=t.start,t.end = t.end,min.t.points = min.t.points,min.distance = min.distance,rank = input$rank=='TRUE',
                                        spline = input$method.intersection=='Spline',spline.df = input$spline.df)

             })

             scores.filtered<-eventReactive(input$filtering,{
               if (is.null(data.exp()) | is.null(mapping()))
                 return()

               t.points.cutoff<-input$t.points.cutoff
               prob.cutoff<-input$prob.cutoff
               dist.cutoff<-input$dist.cutoff
               pval.cutoff<-input$pval.cutoff
               cor.cutoff<-input$cor.cutoff
               x.value.limit<-c(input$x.lower.boundary,input$x.upper.boundary)
               #
               scores.filtered<-score.filter(scores = scores(),t.points.cutoff=t.points.cutoff,prob.cutoff=prob.cutoff,dist.cutoff=dist.cutoff,
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
                                       'left.interval','right.invertal','x.value','y.value','prob','dist','left.pval','right.pval','left.t.points','right.t.points','cor')]


             })


             observeEvent(input$filtering, {
               if (is.null(data.exp()) | is.null(mapping()))
                 return()

               x<- scores.filtered()
               if(nrow(x)>1)
                 x[,-c(1,2,5,6)]<-apply(x[,-c(1,2,5,6)],2,function(x) as.numeric(format(x,digits = 3)))

               rownames(x)<-NULL
               score.show$scores<-x[,c('iso1','iso2','iso1.mean.ratio','iso2.mean.ratio',
                                       'left.interval','right.invertal','x.value','y.value','prob','dist','left.pval','right.pval','left.t.points','right.t.points','cor')]


             })




             output$score.table <- shiny::renderDataTable(options=list(pageLength=10,aoColumnDefs = list(list(sClass="alignLeft",aTargets="_all")) ),{
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
                 write.table(score.show$scores,file,row.names = F,col.names = T,sep=',',quote = F)
               }
             )

             ##plot the density
             # height = 400, width = 600
             output$density <- renderPlotly({
               if(is.null(score.show$scores))
                 return()

               switch.density(x=score.show$scores$x.value,make.plotly = T)
             })


             #########


             g<-eventReactive(input$plot1iso,{
               if(is.null(input$iso1) | is.null(input$iso2) | is.null(data.exp()))
                 return()

               iso1<-input$iso1
               iso2<-input$iso2
               g<-plotTSIS(data2plot = data.exp()[c(iso1,iso2),],
                          iso1 = iso1,
                          iso2 = iso2,
                          scores = scores.filtered(),show.scores = input$show.scores,
                          t.start=input$t.start,
                          t.end=input$t.end,
                          nrep=input$nrep,
                          x.lower.boundary=input$x.lower.boundary,
                          x.upper.boundary=input$x.upper.boundary,
                          prob.cutoff=input$prob.cutoff.switch.points,
                          y.lab = 'Expression',spline=(input$method.intersection=='Spline'),spline.df = input$spline.df,ribbon.plot = (input$ribbon.plot=='Ribbon')

               )
             })

             output$plot.isoform<-renderPlotly({
               if(is.null(g()))
                 return()

               ggplotly(g())
             })
             #

             output$download.1plot <- downloadHandler(

               filename = function() {
                 paste0(input$iso1,' vs ',input$iso2,'.',input$plot1.format)
               },
               content = function(file,format=input$plot1.format) {
                 if(format=='html')
                   suppressWarnings(htmlwidgets::saveWidget(as.widget(ggplotly(g())), file=file,selfcontained=T))
                 else ggsave(file, g(),width = 25,height = 12,units = "cm")
               })


             topn<-eventReactive(input$plotmultiiso,{
               topn<-input$topn
             })

             plotmulti.format<-eventReactive(input$plotmultiiso,{
               plotmulti.format<-input$plotmulti.format
             })



             folderInput <- eventReactive(input$plotmultiiso,{
               x<-paste0(getwd(),'/',input$folder2save)
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




             output$isoform.pairs.to.plot <- shiny::renderDataTable(options=list(pageLength=10,aoColumnDefs = list(list(sClass="alignLeft",aTargets="_all")) ),{
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
                               iso1 = NULL,
                               iso2 = NULL,
                               scores = data2plot(),
                               t.end=input$t.end,
                               t.start = input$t.start,
                               nrep=input$nrep,
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





             #
           }
                                                )
}


