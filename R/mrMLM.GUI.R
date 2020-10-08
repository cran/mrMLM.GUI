mrMLM.GUI<-function(){
  
  server<-(function(input, output,session) {
    
    options(shiny.maxRequestSize=-1)
    
    output$general_ui <- renderUI({
      
      if (input$layout == "mrMLM") {
        
        tabPanel("",
                 
                 h2("Multi-locus GWAS methods",align="center"),
                 br(),
                 column(10,
                        h4("1. Zhang YM, Mao Y, Xie C, Smith H, Luo L, Xu S*. Mapping quantitative trait loci using 
                           naturally occurring genetic variance among commercial inbred lines of maize (", em("Zea mays "),
                           "L.).Genetics 2005;169:2267-2275"),
                        h4("2. Wang SB, Feng JY, Ren WL, Huang B, Zhou L, Wen YJ, Zhang J, Jim M Dunwell, Xu S*, Zhang YM*.
                           Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear
                           model methodology. Scientific Reports 2016;6:19444. (mrMLM)"),
                        h4("3. Tamba CL, Ni YL, Zhang YM*. Iterative sure independence screening EM-Bayesian LASSO algorithm for
                           multi-locus genome-wide association studies. PLoS Computational Biology 2017;13(1):e1005357. (ISIS EM-BLASSO)"),
                        h4("4. Zhang J, Feng JY, Ni YL, Wen YJ, Niu Y, Tamba CL, Yue C, Song QJ, Zhang YM*. pLARmEB: 
                           integration of least angle regression with empirical Bayes for multi-locus genome-wide association
                           studies. Heredity 2017;118(6):517-524. (pLARmEB)"),
                        h4("5. Ren WL, Wen YJ, Jim M Dunwell, Zhang YM*. pKWmEB: integration of Kruskal-Wallis test with
                           empirical Bayes under polygenic background control for multi-locus genome-wide association study. 
                           Heredity 2018;120(3):208-218. (pKWmEB)"),
                        h4("6. Wen YJ, Zhang H, Ni YL, Huang B, Zhang J, Feng JY, Wang SB, Jim M Dunwell, Zhang YM*, Wu R*.
                           Methodological implementation of mixed linear models in multi-locus genome-wide association studies. 
                           Briefings in Bioinformatics 2018;19(4):700-712. (FASTmrEMMA)"),
                        h4("7. Tamba CL, Zhang YM. A fast mrMLM algorithm for multi-locus genome-wide association studies. 
                           bioRxiv, 2018. doi: 10.1101/341784. (FASTmrMLM)"),
                        h4("8. Zhang Yuan-Ming, Jia Zhenyu, Jim M. Dunwell. Editorial: The applications of new multi-locus GWAS 
                            methodologies in the genetic dissection of complex traits. Frontiers in Plant Science 2019, 10: 100."),
                        h4("9. Zhang Ya-Wen, Tamba Cox Lwaka, Wen Yang-Jun, Li Pei, Ren Wen-Long, Ni Yuan-Li, Gao Jun, Zhang Yuan-Ming*. 
                           mrMLM v4.0: An R platform for multi-locus genome-wide association studies. Genomics, Proteomics & Bioinformatics
                           2020; Accept (bioRxiv, 2020.03.04.976464)"),
                        br(),
                        h4("Authors: Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming"),
                        
                        h4("Maintainer: Zhang Yuan-Ming (soyzhang at mail.hzau.edu.cn)"), 
                        
                        h4("mrMLM.GUI version 4.0.2, Realeased October 2020"),
                        
                        offset=1)
                 
                        )
        
      } else if (input$layout == "Start") {
        
        
        fluidPage(  
          titlePanel(""), 
          navlistPanel(
            
            selected="Genotype",
            
            tabPanel("Genotype",
                     
                     fluidRow( 
                       
                       column(12,
                              h3("Genotype"),
                              offset = 1
                       ),
                       
                       column(3,
                              br(),
                              
                              radioButtons("radiomrformat", "Dataset format",
                                           choices = list("mrMLM numeric format", "mrMLM character format","Hapmap (TASSEL) format"),selected ="mrMLM numeric format"),
                              offset = 1
                       ),
                       column(3,
                              br(),
                              fileInput("filegenmr", "Genotypic file",multiple = TRUE,accept = ".csv")
                       ),
                       
                       column(3,
                              br(),
                              radioButtons("dispgenmr", "Display genotype", choices = c(Head = "head",All = "all"),selected = "head"),
                              offset = 1
                       ),
                       column(12,
                              tableOutput("contentsrawgenmr"),
                              offset = 1
                       )
                     ) 
            ) ,
            tabPanel("Phenotype",
                     
                     fluidRow( 
                       column(12,
                              h3("Phenotype"),
                              offset = 1
                       ),
                       column(3,
                              br(),
                              fileInput("filephenmr", "Phenotypic file",multiple = TRUE,accept = ".csv"),
                              offset = 1
                       ),
                       
                       column(3,
                              br(),
                              radioButtons("dispphemr", "Display phenotype", choices = c(Head = "head",All = "all"),selected = "head"),
                              offset = 1
                       ),
                       
                       column(12,
                              br(),
                              tableOutput("contentsrawphemr"),
                              offset = 1
                       )
                     )
            ),   
            
            tabPanel("Kinship",
                     fluidRow(  
                       column(12,
                              h3("Kinship"),
                              offset=1
                       ),
                       column(3,
                              br(),
                              radioButtons("ckmr", "Kinship",
                                           choices = list("Input Kinship (K) matrix file"="inmr", "Calculate Kinship (K) matrix by this software"="cmr"),selected = "cmr"), 
                              conditionalPanel("input.ckmr == 'inmr'",
                                               fileInput("filekmr", "Kinship (K)",multiple = TRUE,accept = ".csv")
                              ),
                              offset = 1
                       ),
                       column(3,
                              br(),
                              radioButtons("dispkinmr", "Display kinship", choices = c(Head = "head",All = "all"),selected = "head"),
                              offset = 1
                       ),
                       
                       column(12,
                              br(),
                              tableOutput("kinshipmr"),
                              
                              offset = 1
                       )
                     )
            ),                 
            
            tabPanel("Population structure",
                     
                     fluidRow(
                       column(12,
                              h3("Population structure"),
                              offset=1
                       ),
                       
                       column(4,
                              br(),
                              radioButtons("instrmr", "Population structure", choices = c("Not included in the model"="notinc" ,"Included"="inc" ),selected = "notinc"), 
                              conditionalPanel("input.instrmr == 'inc'",
                                               radioButtons("poptypeall", "Population structure type", choices = c("Q matrix"="moQ" ,"Main principal components"="pcaQ","Evolutionary population structure"="evoQ"),selected = "moQ")),
                              offset = 1
                       ),
                       column(4,
                              br(),
                              radioButtons("dispstrmr", "Display population structure", choices = c(Head = "head",All = "all"),selected = "head"),
                              conditionalPanel("input.instrmr == 'inc'",
                                               fileInput("filestrmr", "Population structure",multiple = TRUE,accept = ".csv")),
                              offset = 1
                       ),
                       
                       column(12,
                              br(),
                              tableOutput("psmr"),
                              offset = 1
                       )
                     )
                     
            ),
            
            
            tabPanel("Covariate",
                     fluidRow(
                       column(12,
                              h3("Covariate"),
                              offset=1
                       ),
                       column(4,
                              br(),
                              radioButtons("incov", "Covariate", choices = c("Not included in the model"="nocov" ,"Included"="cov" ),selected = "nocov"), 
                              conditionalPanel("input.incov == 'cov'",
                                               fileInput("filecov", "Covariate",multiple = TRUE,accept = ".csv")),
                              
                              offset = 1
                       ),
                       column(4,
                              br(),
                              radioButtons("dispcov", "Display covariate", choices = c(Head = "head",All = "all"),selected = "head"),
                              offset = 1
                       ),
                       column(12,
                              br(),
                              tableOutput("covcontent"),
                              offset = 1
                       )
                     )
                     
            ),
            
            tabPanel("Method select & Parameter settings",
                     
                     fluidRow(  
                       column(4,
                              br(),
                              checkboxGroupInput("Method", "Method selection",
                                                 c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),selected=c("mrMLM","FASTmrMLM")),
                              conditionalPanel("input.Method.indexOf('mrMLM') > -1||input.Method.indexOf('FASTmrMLM') > -1",
                                               textInput("scgmr", label = "Search radius (kb) of candidate gene (mrMLM & FASTmrMLM):",value = "20")),
                              
                              conditionalPanel("input.Method.indexOf('FASTmrEMMA') > -1",
                                               radioButtons("chdofunFME", "Likelihood Function (FASTmrEMMA)", choices = c("REML","ML"))),
                              conditionalPanel("input.Method.indexOf('pLARmEB') > -1",
                                               textInput("scgPLA", label = "No. of potentially associated variables selected by LARS (pLARmEB):",value = "50"),
                                               radioButtons("BootPLA", "Bootstrap (pLARmEB)", choices = c("TRUE","FALSE"),selected="FALSE")
                              ),
                              
                              
                              offset=1 
                       ),
                       column(4,
                              br(),
                              textInput("lodmr", label = "Critical LOD score (All methods)",value = "3"),
                              textInput("SavePath", "Save path",value = "C:/Users/Administrator/Desktop"),
                              textInput("trait", "Traits analyzed", value="1"),
                              radioButtons("Drawall", "Draw plot (All methods)", choices = c("TRUE","FALSE"),selected="FALSE"),
                              conditionalPanel("input.Drawall=='TRUE'",
                                               radioButtons("Plotformat", "Plot format (All methods)", choices = c("*.png", "*.tiff", "*.jpeg","*.pdf"))
                              )
                              
                       ),
                       column(12, 
                              br(),
                              actionButton("runmr", label = "Run",icon("paper-plane"), 
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4",width=280),
                              offset = 1
                       ),
                       column(12,
                              uiOutput("result")
                       )
                       
                     )
            ),
            
            
            tabPanel("Manhattan Plot",id="manp",
                     
                     fluidRow(
                       
                       column(6,
                              fileInput("filemanIn", "Intermediate result to draw Manhattan plot",multiple = TRUE,accept = ".csv"),
                              fileInput("filemanFin", "Final result to draw Manhattan plot",multiple = TRUE,accept = ".csv"),
                              radioButtons("plmanformat", "Manhattan plot format", c("*.png", "*.tiff", "*.jpeg","*.pdf")),
                              conditionalPanel("input.plmanformat != '*.pdf'",
                                               textInput("manwidthmr", label = "Figure width (px):",value = "28000"),
                                               textInput("manheimr", label = "Figure height (px)",value = "7000"),
                                               textInput("manwordremr", label = "Word resolution (1/72 inch, ppi):",value = "60"),
                                               textInput("manfigureremr", label = "Figure resolution (ppi):",value = "600")
                              ),
                              conditionalPanel("input.plmanformat == '*.pdf'",
                                               textInput("manwidthpdf", label = "Figure width (inches):",value = "16"),
                                               textInput("manheipdf", label = "Figure height (inches)",value = "4"),
                                               textInput("manpoipdf", label = "Word resolution (1/72 inch, ppi):",value = "20")                
                              ),
                              
                              textInput("labelsize", label = "Size of all the three labels",value = "0.8"),
                              textInput("TckLwd", label = "Size of scale values",value = "0.7"),
                              textInput("CoorLwd", label = "Width of all the three axes",value = "5"),
                              textInput("TckLen", label = "Length of tick marks",value = "-0.03")
                           ),
                       column(6,
                             textInput("VerLabDis", label = "Distance between label and axis",value = "1.5"),
                              #textInput("HorLabDis", label = "Distance between label and horizontal axis",value = "1.5"),
                              textInput("VerTckDis", label = "Distance between scale values and tick marks",value = "0.4"),
                              #textInput("HorTckDis", label = "Distance between scale values and horizontal tick marks",value = "0.2"),
                              textInput("logtimes", label = "Magnification of {-log10(P-value)}",value = "2"),
                              textInput("LODtimes", label = "Magnification of {LOD score}",value = "1.2"),
                              
                              radioButtons("markgene", "Mark genes or not", c("TRUE", "FALSE"),selected = "FALSE"),
                              conditionalPanel("input.markgene == 'TRUE'",
                                               textInput("posX", label = "Numeric vectors of x axis",value = "139,195"),
                                               textInput("posY", label = "Numeric vectors of y axis",value = "7.5,7"),
                                               textInput("Genename",label = "Character vectors of gene names that mark in the plot ",value = "Gene1,Gene2"),
                                               textInput("colorname",label = "Colour of gene names",value = "blue")                
                                               
                              ),
                              
                              textInput("Saveplot", "Save path",value = "C:/Users/Administrator/Desktop"),
                              br(),

                              actionButton("drawman", label = "Draw Manhattan plot",icon("paper-plane"), 
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4",width="55%"),
                             plotOutput("manout")
                         )
                      )
            ),
            
            tabPanel("QQ Plot",id="qqp",
                     
                     fluidRow(
                       
                       column(6,
                              fileInput("fileqqmr", "File to draw QQ plot",multiple = TRUE,accept = ".csv"),
                              radioButtons("plqqformat", "QQ plot format", c("*.png", "*.tiff", "*.jpeg","*.pdf")),
                              conditionalPanel("input.plqqformat != '*.pdf'",
                                               textInput("qqwidthmr", label = "Figure width (px):",value = "10000"),
                                               textInput("qqheimr", label = "Figure height (px)",value = "10000"),
                                               textInput("qqwordremr", label = "Word resolution (1/72 inch, ppi):",value = "60"),
                                               textInput("qqfigureremr", label = "Figure resolution (ppi):",value = "600")
                              ),
                              conditionalPanel("input.plqqformat == '*.pdf'",
                                               textInput("qqwidthpdf", label = "Figure width (inches):",value = "7"),
                                               textInput("qqheipdf", label = "Figure height (inches)",value = "7"),
                                               textInput("qqpoipdf", label = "Word resolution (1/72 inch, ppi):",value = "25")                
                              ),
                              textInput("labelsizeqq", label = "Size of all the two labels",value = "0.7"),
                              textInput("CoorLwdqq", label = "Size of all the two axes",value = "3")
                              
                       ),
                       
                       column(6,
                              textInput("TckLenqq", label = "Length of tick marks",value = "-0.02"),
                              textInput("TckLwdqq", label = "Size of scale values",value = "0.6"),
                              textInput("VerLabDisqq", label = "Distance between label and vertical axis",value = "1.1"),
                              textInput("HorLabDisqq", label = "Distance between label and horizontal axis",value = "1"),
                              textInput("VerTckDisqq", label = "Distance between scale values and vertical tick marks",value = "0.3"),
                              textInput("HorTckDisqq", label = "Distance between scale values and horizontal tick marks",value = "0.02"),
                              textInput("pstand", label = "Critical P-value of deleting points",value = "0.9"),
                              textInput("Saveplotqq", "Save path",value = "C:/Users/Administrator/Desktop"),
                              br(),
                              br(),
                              actionButton("drawqq", label = "Draw QQ plot",icon("paper-plane"), 
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4",width="55%"),
                              
                              plotOutput("qqout")
                              
                       )
                     ) 
            )
          ),
          
          br(),
          actionButton("manl", label = "User manual",width=280,
                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
          uiOutput("manl")
          
        )
      }
      
    })
    
    
    Regen<-reactive({
      req(input$metagen)
      datagenRaw<-as.matrix(fread(input$metagen$datapath,header = FALSE,stringsAsFactors=T))
    })
    
    
    Rephe<-reactive({
      req(input$metaphe)
      datapheRaw<-as.matrix(fread(input$metaphe$datapath,header = FALSE,stringsAsFactors=T))
    })
    
    
    Repop<-reactive({
      req(input$filemetastr)
      dataQRaw<-as.matrix(fread(input$filemetastr$datapath,header = FALSE,stringsAsFactors=T))
    })
    
    Recov<-reactive({
      req(input$filemetacov)
      datacovRaw<-as.matrix(fread(input$filemetacov$datapath,header = FALSE,stringsAsFactors=T))
    })
    
    
    RedatamrMLM<-reactive({
      req(input$metafilemrMLM)
      datamrMLM<-fread(input$metafilemrMLM$datapath)
    })
    
    Redataother<-reactive({
      req(input$metafileother)
      dataother<-fread(input$metafileother$datapath)
    })
    
    
    
    ReadDataMeta<-reactive({
      psmatrixRaw<-NULL
      covmatrixRaw<-NULL
      genRaw<-Regen()
      pheRaw1q<-Rephe()
      if(input$metainps=="metastr"&&input$filemetastr !=""){
        psmatrixRaw<-Repop()
      }
      if(input$metacov=="metahavecov"&&input$filemetacov !=""){
        covmatrixRaw<-Recov()
      }
      phename<-as.matrix(pheRaw1q[1,2:ncol(pheRaw1q)])
      output<-list(genRaw=genRaw,pheRaw1q=pheRaw1q,psmatrixRaw=psmatrixRaw,covmatrixRaw=covmatrixRaw,phename=phename)
    })
    
    
    MetaAna<-reactive({
      
      DoData<-function(genRaw=NULL,Genformat=NULL,pheRaw1q=NULL,kkRaw=NULL,psmatrixRaw=NULL,covmatrixRaw=NULL,trait=NULL,type=NULL,PopStrType=NULL){
        inputform<-Genformat
        pheRaw1qq<-as.matrix(pheRaw1q[,2:ncol(pheRaw1q)])
        pheRaw1<-cbind(pheRaw1q[,1],pheRaw1qq[,trait])
        pheRaw2<-pheRaw1[-1,]
        pheRaw3<-as.data.frame(pheRaw2,stringsAsFactors=FALSE)
        pheRaw4<-as.matrix(pheRaw3[is.na(pheRaw3[,2])==F,])
        pheRawthem<-matrix(c(pheRaw1[1,1]," "),1,)
        pheRaw<-rbind(pheRawthem,pheRaw4)
        row.names(pheRaw)<-NULL
        pheRaw<-as.matrix(pheRaw)
        
        if(type==1&&inputform==1){
          genRawz<-genRaw[-1,-c(1:4)]
          genRawz2<-gsub("0","0.5",genRawz)
          genRawz3<-gsub("-1","0",genRawz2)
          genRawz4<-cbind(genRaw[-1,c(1:4)],genRawz3)
          genRaw<-rbind(genRaw[1,],genRawz4)
        }else{
          genRaw<-genRaw
        }
        if(inputform==1){
          nameGen <- as.matrix(genRaw[1,],1,)
          namePhe <- as.matrix(pheRaw[,1],,1)
          sameName <- intersect(nameGen,namePhe)
          ##########To find the location of the same name
          locGen <- match(sameName,nameGen)
          locPhe <- match(sameName,namePhe)
          ##########Produce new genotype matrix and phenotype matrix
          hapName <- matrix(c("rs#","chrom","pos","genotype for code 1"),1,)
          hapHave <- intersect(nameGen,hapName)
          locHap <- match(hapHave,nameGen)
          newGenloc <- c(locHap,locGen)
          newPheloc <- locPhe
          newGen <- as.matrix(genRaw[-1,newGenloc])
          newPhe <- as.matrix(pheRaw[newPheloc,])
          nnhap <- length(hapHave)
          rownewGen <- dim(newGen)[1]
          colnewGen <- dim(newGen)[2]
          rownewPhe <- dim(newPhe)[1]
          ###########To show on the table ----newGen
          newGen <-rbind(genRaw[1,newGenloc],newGen)
          ###########To be computed ----gen
          locChr <- as.numeric(which(newGen[1,]=="chrom"))
          locPos <- as.numeric(which(newGen[1,]=="pos"))
          needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
          needGen <- newGen[,needloc]
          
          genq<-as.matrix(needGen[-1,])
          gen<-big.matrix(nrow(genq),ncol(genq),type="double",shared = FALSE)
          gen[,]<-genq[,]
          rm(newGen,needGen,genq)
          gc()
          ###########To show on the table ----newPhe
          pheRaw[1,2]<-"  "
          newPhe<-rbind(pheRaw[1,],newPhe)
          ###########To be computed ----phe
          phe<-as.matrix(newPhe[-1,-1])
          phe<-matrix(as.numeric(phe),nrow=nrow(phe))
          outATCG<-NULL
        }else if(inputform==2){
          ##########To find the same individual ID between genotype and phenotype
          nameGen <- as.matrix(genRaw[1,],1,)
          namePhe <- as.matrix(pheRaw[,1],,1)
          sameName <- intersect(nameGen,namePhe)
          ##########To find the location of the same name
          locGen <- match(sameName,nameGen)
          locPhe <- match(sameName,namePhe)
          ##########Produce new genotype matrix and phenotype matrix
          hapName <- matrix(c("rs#","chrom","pos"),1,)
          hapHave <- intersect(nameGen,hapName)
          locHap <- match(hapHave,nameGen)
          newGenloc <- c(locHap,locGen)
          newPheloc <- locPhe
          newGen <- as.matrix(genRaw[-1,newGenloc])
          newPhe <- as.matrix(pheRaw[newPheloc,])
          ##########Transfer ATCG to numeric
          nnhap <- length(hapHave)
          rownewGen <- dim(newGen)[1]
          colnewGen <- dim(newGen)[2]
          rownewPhe <- dim(newPhe)[1]
          computeGen <- newGen[,(nnhap+1):colnewGen]
          colComGen <- ncol(computeGen)
          referSam <- as.vector(computeGen[,1])
          ATCGloc <- c(which(computeGen[,1]=="A"),which(computeGen[,1]=="T"),which(computeGen[,1]=="C"),which(computeGen[,1]=="G"))
          NNRRloc <- setdiff(c(1:rownewGen),ATCGloc)
          for(i in 2:colComGen)
          {
            if(length(NNRRloc)>0){
              referSam[NNRRloc] <- as.vector(computeGen[NNRRloc,i])
              ATCGlocLoop <- c(which(computeGen[NNRRloc,i]=="A"),which(computeGen[NNRRloc,i]=="T"),which(computeGen[NNRRloc,i]=="C"),which(computeGen[NNRRloc,i]=="G"))
              NNRRloc <- setdiff(NNRRloc,NNRRloc[ATCGlocLoop])
            }else{
              break
            }
          }
          for(i in 1:rownewGen)
          {
            tempSel1 <- as.vector(c(which(computeGen[i,]=="A"),which(computeGen[i,]=="T"),which(computeGen[i,]=="C"),which(computeGen[i,]=="G")))
            tempSel2 <- as.vector(c(which(computeGen[i,]==referSam[i])))
            notRef <- setdiff(tempSel1,tempSel2)
            notATCG <- setdiff(c(1:colComGen),tempSel1)
            computeGen[i,tempSel2] <- as.numeric(1)
            
            if(type==1){
              computeGen[i,notRef] <- as.numeric(0)
              computeGen[i,notATCG] <- as.numeric(0.5)
            }else{
              computeGen[i,notRef] <- as.numeric(-1)
              computeGen[i,notATCG] <- as.numeric(0)
            }
          }
          outATCG<-as.matrix(referSam)
          ###########To show on the table ----newGen
          newGen <- cbind(newGen[,1:nnhap],computeGen)
          newGen <-rbind(genRaw[1,newGenloc],newGen)
          rm(computeGen)
          gc()
          ###########To be computed ----gen
          locChr <- as.numeric(which(newGen[1,]=="chrom"))
          locPos <- as.numeric(which(newGen[1,]=="pos"))
          needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
          needGen<-newGen[,needloc]
          
          genq<-as.matrix(needGen[-1,])
          gen<-big.matrix(nrow(genq),ncol(genq),type="double",shared = FALSE)
          gen[,]<-genq[,]
          rm(newGen,needGen,genq)
          gc()
          ###########To show on the table ----newPhe
          pheRaw[1,2]<-"  "
          newPhe<-rbind(pheRaw[1,],newPhe)
          ###########To be computed ----phe
          phe<-as.matrix(newPhe[-1,-1])
          phe<-matrix(as.numeric(phe),nrow=nrow(phe))
        }else if(inputform==3){
          ##########To find the same individual ID between genotype and phenotype
          nameGen<-as.matrix(genRaw[1,],1,)
          namePhe<-as.matrix(pheRaw[,1],,1)
          sameName<-intersect(nameGen,namePhe)
          ##########To find the location of the same name
          locGen<-match(sameName,nameGen)
          locPhe<-match(sameName,namePhe)
          ##########Produce new genotype matrix and phenotype matrix
          hapName<-matrix(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode"),1,)
          hapHave<-intersect(nameGen,hapName)
          locHap<-match(hapHave,nameGen)
          newGenloc<-c(locHap,locGen)
          newPheloc<-locPhe
          newGen<-as.matrix(genRaw[-1,newGenloc])
          newPhe<-as.matrix(pheRaw[newPheloc,])
          ##########Transfer ATCG to numeric
          nnhap<-length(hapHave)
          rownewGen<-dim(newGen)[1]
          colnewGen<-dim(newGen)[2]
          rownewPhe<-dim(newPhe)[1]
          computeGen<-newGen[,(nnhap+1):colnewGen]
          colComGen<-ncol(computeGen)
          referSam<-as.vector(computeGen[,1])
          ATCGloc<-c(which(computeGen[,1]=="AA"),which(computeGen[,1]=="TT"),which(computeGen[,1]=="CC"),which(computeGen[,1]=="GG"))
          NNRRloc<-setdiff(c(1:rownewGen),ATCGloc)
          for(i in 2:colComGen)
          {
            if(length(NNRRloc)>0){
              referSam[NNRRloc]<-as.vector(computeGen[NNRRloc,i])
              ATCGlocLoop<-c(which(computeGen[NNRRloc,i]=="AA"),which(computeGen[NNRRloc,i]=="TT"),which(computeGen[NNRRloc,i]=="CC"),which(computeGen[NNRRloc,i]=="GG"))
              NNRRloc<-setdiff(NNRRloc,NNRRloc[ATCGlocLoop])
            }else{
              break
            }
          }
          for(i in 1:rownewGen)
          {
            tempSel1<-as.vector(c(which(computeGen[i,]=="AA"),which(computeGen[i,]=="TT"),which(computeGen[i,]=="CC"),which(computeGen[i,]=="GG")))
            tempSel2<-as.vector(c(which(computeGen[i,]==referSam[i])))
            notRef<-setdiff(tempSel1,tempSel2)
            notATCG<-setdiff(c(1:colComGen),tempSel1)
            computeGen[i,tempSel2]<-as.numeric(1)
            if(type==1){
              computeGen[i,notRef]<-as.numeric(0)
              computeGen[i,notATCG]<-as.numeric(0.5)
            }else{
              computeGen[i,notRef]<-as.numeric(-1)
              computeGen[i,notATCG]<-as.numeric(0)
            }
          }
          outATCG<-as.matrix(referSam)
          ###########To show on the table ----newGen
          newGen<-cbind(newGen[,1:nnhap],computeGen)
          newGen<-rbind(genRaw[1,newGenloc],newGen)
          
          rm(computeGen)
          gc()
          ###########To be computed ----gen
          locChr<-as.numeric(which(newGen[1,]=="chrom"))
          locPos<-as.numeric(which(newGen[1,]=="pos"))
          needloc<-c(locChr,locPos,(nnhap+1):colnewGen)
          needGen<-newGen[,needloc]
          
          genq<-as.matrix(needGen[-1,])
          gen<-big.matrix(nrow(genq),ncol(genq),type="double",shared = FALSE)
          gen[,]<-genq[,]
          rm(newGen,needGen,genq)
          gc()
          ###########To show on the table ----newPhe
          pheRaw[1,2]<-"  "
          newPhe<-rbind(pheRaw[1,],newPhe)
          ###########To be computed ----phe
          phe<-as.matrix(newPhe[-1,-1])
          phe<-matrix(as.numeric(phe),nrow=nrow(phe))
          
        }
        
        if(is.null(kkRaw)){
          kk<-NULL
        }else{
          # kkPre<-as.matrix(kkRaw[-1,-1])
          # nameKin<-as.matrix(kkRaw[-1,1])
          kkPre<-as.matrix(kkRaw[,-1])
          nameKin<-as.matrix(kkRaw[,1])
          sameGenKin<-intersect(sameName,nameKin)
          locKin<-match(sameGenKin,nameKin)
          kk<-kkPre[locKin,locKin]
          kk<-matrix(as.numeric(kk),nrow=nrow(kk))
        }
        
        if(is.null(psmatrixRaw)){
          psmatrix<-NULL
        }else{
          nnpprow<-dim(psmatrixRaw)[1]
          nnppcol<-dim(psmatrixRaw)[2]
          psmatrixRaw[1,2:nnppcol]<-"  "
          psmatrixPre<-psmatrixRaw[3:nnpprow,]
          namePop<-as.matrix(psmatrixPre[,1])
          sameGenPop<-intersect(sameName,namePop)
          locPop<-match(sameGenPop,namePop)
          selectpsmatrixq<-psmatrixPre[locPop,-1]
          
          if(PopStrType=="metamoQ"||PopStrType=="moQ"){
            selectpsmatrix<-matrix(as.numeric(selectpsmatrixq),nrow = length(locPop))
            coldelet<-which.min(apply(selectpsmatrix,2,sum))
            psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
          }else if(PopStrType=="metapcaQ"||PopStrType=="pcaQ"){
            psmatrix<-matrix(as.numeric(selectpsmatrixq),nrow = length(locPop))
          }else if(PopStrType=="metaevoQ"||PopStrType=="evoQ"){
            otrait_ind<-sort(unique(selectpsmatrixq))
            pop_col<-length(otrait_ind)-1
            pop_each<-numeric()
            for(j in 1:length(selectpsmatrixq)){
              if(selectpsmatrixq[j]==otrait_ind[1]){
                pop_0<-matrix(-1,1,pop_col)
              }else{
                pop_0<-matrix(0,1,pop_col)
                popnum_loc<-which(otrait_ind[]==selectpsmatrixq[j])
                pop_0[1,popnum_loc-1]<-1
              }
              pop_each<-rbind(pop_each,pop_0)
            }
            psmatrix=pop_each
          }
          
        }
        if(is.null(covmatrixRaw)){
          phe<-phe
        }else{
          nncovrow<-nrow(covmatrixRaw)
          covmatrixPre<-covmatrixRaw[3:nncovrow,]
          namecov<-as.matrix(covmatrixPre[,1])
          sameGencov<-intersect(sameName,namecov)
          loccov<-match(sameGencov,namecov)
          selectcovmatrixq<-covmatrixPre[loccov,-1]
          covname<-covmatrixRaw[2,-1]
          label<-substr(covname,1,3)
          if(("Cat"%in%label)&&("Con"%in%label)){
            cat_loc<-as.numeric(which(label=="Cat"))
            con_loc<-as.numeric(which(label=="Con"))
            selectcovmatrixqq<-selectcovmatrixq
            selectcovmatrixq<-selectcovmatrixq[,cat_loc]
            covnum<-t(selectcovmatrixq)
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){
              otrait_ind<-sort(unique(covnum[i,]))
              cov_col<-length(otrait_ind)-1
              col_each<-numeric()
              for(j in 1:length(covnum[i,])){
                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)
                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)
                
              }
              yygg1<-cbind(yygg1,col_each)
            }
            yygg1<-cbind(yygg1,as.matrix(selectcovmatrixqq[,con_loc]))
          }else if(all(label=="Cat")){
            covnum<-t(selectcovmatrixq)
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){
              otrait_ind<-sort(unique(covnum[i,]))
              cov_col<-length(otrait_ind)-1
              col_each<-numeric()
              for(j in 1:length(covnum[i,])){
                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)
                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)
                
              }
              yygg1<-cbind(yygg1,col_each)
            }
          }else if(all(label=="Con")){
            yygg1<-selectcovmatrixq
          }
          
          W.orig<-matrix(1,nrow(phe),1)
          xenvir<-cbind(W.orig,yygg1)
          xenvir<-apply(xenvir,2,as.numeric)
          beta<-solve(t(xenvir)%*%xenvir)%*%t(xenvir)%*%phe
          phe<-phe-xenvir%*%beta+W.orig
        }
        genRaw<-genRaw[,1:12]
        doresult<-list(gen=gen,phe=phe,outATCG=outATCG,genRaw=genRaw,kk=kk,psmatrix=psmatrix)
        
        return(doresult)
        
      }
      
      inputData<-function(readraw,Genformat=NULL,method=NULL,trait=NULL,PopStrType=NULL){
        
        doMR<-NULL;doFME<-NULL
        
        if("mrMLM"%in%method){
          doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
        }
        
        if("FASTmrMLM"%in%method){
          doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType)
        }  
        
        if("FASTmrEMMA"%in%method){
          doFME<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=1,PopStrType)  
        }
        
        if("pLARmEB"%in%method){
          doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType) 
        }
        if("pKWmEB"%in%method){
          doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType) 
        }  
        
        if("ISIS EM-BLASSO"%in%method){
          doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,readraw$covmatrixRaw,trait,type=2,PopStrType) 
        }
        
        output<-list(doMR=doMR,doFME=doFME) 
        return(output)
        
      }
      
      output=list(Dodata=DoData,inputData=inputData)
    })
    
    
    
    
    
    manual<-eventReactive(input$manl,{
      
      RShowDoc("Instruction",package="mrMLM.GUI") 
      
    })
    
    output$manl<-renderUI(manual())
    
    
    genRawmr<-reactive({
      req(input$filegenmr)
      grm<-fread(input$filegenmr$datapath,header = FALSE,stringsAsFactors=T)
      grm<-as.matrix(grm)
    })
    
    CoRaGenmr<-reactive({
      genRaw<-genRawmr()
      showgenRaw<-genRaw[-1,]
      colnames(showgenRaw)<-genRaw[1,]
      showgenRaw<-as.data.frame(showgenRaw)
    })
    
    
    output$contentsrawgenmr <- renderTable({
      if(input$dispgenmr == "head"){
        return(head(CoRaGenmr()))
      }else{
        return(CoRaGenmr())
      }
    })
    
    
    pheRawmr1<-reactive({
      req(input$filephenmr)
      prm1<-as.matrix(fread(input$filephenmr$datapath,header=FALSE,stringsAsFactors=T)) 
    }) 
    
    CoRaPhemr<-reactive({
      pheRaw1<-pheRawmr1()
      pheRaw2<-pheRaw1[-1,]
      pheRaw3<-as.data.frame(pheRaw2)
      colnames(pheRaw3)<-pheRaw1[1,]
      pheRaw3
    })
    
    output$contentsrawphemr <- renderTable({
      
      if(input$dispphemr == "head"){
        return(head(CoRaPhemr()))
      }else{
        return(CoRaPhemr())
      } 
    })
    
    
    kinmrInput<-reactive({
      req(input$filekmr)
      kkRaw<-as.matrix(fread(input$filekmr$datapath,header = FALSE,stringsAsFactors=T))
    })
    
    
    kinmrshow<-reactive({
      kkRaw<-kinmrInput()
      nnkk<-dim(kkRaw)[1]
      kkRaw[1,2:nnkk]<-" "
      kkRaw2<-kkRaw[-1,]
      kkRaw3<-as.data.frame(kkRaw2)
      colnames(kkRaw3)<-kkRaw[1,]
      kkRaw3
      
    })
    
    output$kinshipmr<-renderTable({
      
      if(input$ckmr=="inmr"){
        
        if(input$dispkinmr == "head"){
          return(head(kinmrshow()))
        }else{
          return(kinmrshow())
        }
        
      }
    })
    
    observe({
      shinyjs::toggleState("filekmr",input$ckmr!="cmr")
    })
    
    
    strcmr<-reactive({
      req(input$filestrmr)
      psmatrixRaw<-as.matrix(fread(input$filestrmr$datapath,header = FALSE,stringsAsFactors=T))
      
    })
    
    
    covmr<-reactive({
      req(input$filecov)
      covmatrixRaw<-as.matrix(fread(input$filecov$datapath,header = FALSE,stringsAsFactors=T))
    })
    
    
    strcmrshow<-reactive({
      psmatrixRaw<-strcmr()
      nnppcol<-dim(psmatrixRaw)[2]
      psmatrixRaw[1,2:nnppcol]<-"  "
      psmatrixRaw2<-psmatrixRaw[-1,]
      psmatrixRaw3<-as.data.frame(psmatrixRaw2)
      colnames(psmatrixRaw3)<-psmatrixRaw[1,]
      psmatrixRaw3
    })
    
    
    output$psmr<-renderTable({
      if(input$instrmr=="inc"){
        if(input$dispstrmr == "head"){
          return(head(strcmrshow()))
        }else{
          return(strcmrshow())
        }
        
      } 
      
    })
    
    covmrshow<-reactive({
      covmatrixRaw<-covmr()
      nncovcol<-dim(covmatrixRaw)[2]
      covmatrixRaw[1,2:nncovcol]<-"  "
      covmatrixRaw2<-covmatrixRaw[-1,]
      covmatrixRaw3<-as.data.frame(covmatrixRaw2)
      colnames(covmatrixRaw3)<-covmatrixRaw[1,]
      covmatrixRaw3
    })
    
    output$covcontent<-renderTable({
      if(input$incov=="cov"){
        if(input$dispcov == "head"){
          return(head(covmrshow()))
        }else{
          return(covmrshow())
        }
        
      } 
      
    })
    
    observe({
      shinyjs::toggleState("filestrmr",input$instrmr!="notinc")
    })
    
    
    ReadData<-reactive({
      kkRaw<-NULL
      psmatrixRaw<-NULL
      covmatrixRaw<-NULL
      
      genRaw<-genRawmr()
      pheRaw1q<- pheRawmr1()
      
      if(input$ckmr=="inmr"&&input$filekmr!=""){
        kkRaw<-kinmrshow()
      }
      if(input$instrmr=="inc"&&input$filestrmr !=""){
        psmatrixRaw<-strcmr()
      }
      
      if(input$incov=="cov"&&input$filecov !=""){
        covmatrixRaw<-covmr()
      }
      phename<-as.matrix(pheRaw1q[1,2:ncol(pheRaw1q)])
      output<-list(genRaw=genRaw,pheRaw1q=pheRaw1q,kkRaw=kkRaw,psmatrixRaw=psmatrixRaw,covmatrixRaw=covmatrixRaw,phename=phename)
    })
    
    
    MR<-eventReactive(input$runmr,{
      dir<-input$SavePath
      metainput<-MetaAna()
      DoData<-metainput[[1]]
      InputDataF<-metainput[[2]]
      #meta_analysis<-metainput[[4]]
      svrad<-as.numeric(input$scgmr)
      svmlod<-as.numeric(input$lodmr)
      lars1<-as.numeric(input$scgPLA)
      Likelihood<-input$chdofunFME
      Bootstrap<-input$BootPLA 
      DrawPlot<-input$Drawall
      Plotformat<-input$Plotformat
      Resolution<-input$Resolution
      PopStrType<-input$poptypeall
      
      if(DrawPlot==TRUE){
        manhattan_mrMLM<-function(data_in,data_fin,mar=c(2.9,2.8,0.7,2.8),VerLabDis=1.5,HorLabDis=1.5,
                                  HorTckDis=0.2,VerTckDis=0.4,label_size=0.8,CoorLwd=5,
                                  TckLen=-0.03,TckLwd=0.7,log_times=2,LOD_times=1.2,lodline){
          
          ###########Data process#################
          ###########intermediate result
          method<-unique(data_in[,3])
          data_method<-list(NULL)
          for(i in 1:length(method)){
            data_method[[i]]<-data_in[which(data_in[,3]==method[i]),]
          }
          logp_4method<-numeric()
          for(i in 1:length(method)){
            method_p<-data_method[[i]][,8]
            logp_4method<-cbind(logp_4method,method_p) 
          }
          logp_4method<-apply(logp_4method,2,as.numeric)
          p_4method<-10^-logp_4method
          p_median<-apply(p_4method,1,median)
          locsub<-which(p_median==0)
          pmin<-min(p_median[p_median!=0])
          subvalue<-10^(1.1*log10(pmin))
          p_median[locsub]<-subvalue
          data_p<-as.matrix(p_median)
          data_num<-as.matrix(seq(1:length(p_median)))
          data_chr<-as.matrix(data_method[[1]][,5])
          data_pos<-as.matrix(data_method[[1]][,6])
          manresult<-cbind(data_chr,data_pos,data_p,data_num)
          manresult<-apply(manresult,2,as.numeric)
          colnames(manresult)<-c("Chromosome","BPnumber","P-value","SNPname")
          manresult<-as.data.frame(manresult)
          #######final result##################
          data_fin_method<-unique(data_fin[,3])
          data_fin_method_length<-1:length(unique(data_fin[,3]))
          for(r in 1:length(unique(data_fin[,3]))){
            data_fin[which(data_fin[,3]==data_fin_method[r]),3]<-r
          }
          data_fin_mark<-matrix(data_fin[,c(5,6,8,3)],,4)
          data_fin_mark<-matrix(apply(data_fin_mark,2,as.numeric),,4)
          data_fin_mark_chr<-matrix(data_fin_mark[order(data_fin_mark[,1]),],,4)
          data_fin_mark_order<-numeric()
          for(i in c(unique(data_fin_mark_chr[,1]))){
            data_fin_mark_erery_chr<-matrix(data_fin_mark_chr[which(data_fin_mark_chr[,1]==i),],,4)
            data_fin_mark_pos<-matrix(data_fin_mark_erery_chr[order(data_fin_mark_erery_chr[,2]),],,4)
            all_pos<-unique(data_fin_mark_pos[,2])
            all_pos_maxlod<-numeric()
            for(ii in 1:length(all_pos)){
              all_pos_every<-matrix(data_fin_mark_pos[which(data_fin_mark_pos[,2]==all_pos[ii]),],,4)
              lod_me<-median(all_pos_every[,3])
              all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,all_pos_every[1,4])
              if(nrow(all_pos_every)>=2){
                all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,max(data_fin_mark[,4])+1)
              }
              all_pos_maxlod<-rbind(all_pos_maxlod,all_pos_every_median)
            }
            data_fin_mark_order<-rbind(data_fin_mark_order,all_pos_maxlod)
          }
          snpOfInterest<-numeric()
          for(i in c(unique(data_fin_mark_order[,1]))){
            manresult_chr<-manresult[which(manresult[,1]==i),]
            data_fin_mark_order_chr<-matrix(data_fin_mark_order[which(data_fin_mark_order[,1]==i),],,4)
            mark_loc<-manresult_chr[which(manresult_chr[,2]%in%data_fin_mark_order_chr[,2]),4]
            snpOfInterest<-c(snpOfInterest,mark_loc) 
          }
          bpnumber <- numeric()
          chrnum <- unique(manresult[,1])
          for(i in 1:length(chrnum))
          {
            bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(manresult[,1]==chrnum[i])))))
          }
          manresult2<-cbind(manresult[,1],bpnumber,manresult[,3:4])
          colnames(manresult2)<-c("Chromosome","BPnumber","P-value","SNPname")
          ##########prepare for data#############################
          x<-manresult2;col=c("lightgreen","lightskyblue");logp=TRUE
          chr = "Chromosome";bp ="BPnumber";p ="P-value";snp="SNPname";
          highlight<-snpOfInterest
          CHR=BP=P=index=NULL
          d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
          if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
          d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
          d <- d[order(d$CHR, d$BP), ]
          if (logp) {
            d$logp <- -log10(d$P)
          } else {
            d$logp <- d$P
          }
          d$pos=NA
          d$index=NA
          ind = 0
          for (i in unique(d$CHR)){
            ind = ind + 1
            d[d$CHR==i,]$index = ind
          }
          
          nchr = length(unique(d$CHR))
          if (nchr==1) { ## For a single chromosome
            ## Uncomment the next two linex to plot single chr results in Mb
            #options(scipen=999)
            #d$pos=d$BP/1e6
            d$pos=d$BP
            ticks=floor(length(d$pos))/2+1
            xlabel = paste('Chromosome',unique(d$CHR),'position')
            labs = ticks
          } else { ## For multiple chromosomes
            lastbase=0
            ticks=NULL
            for (i in unique(d$index)) {
              if (i==1) {
                d[d$index==i, ]$pos=d[d$index==i, ]$BP
              } else {
                lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
                d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
              }
              # Old way: assumes SNPs evenly distributed
              # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
              # New way: doesn't make that assumption
              ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
            }
            xlabel = 'Chromosomes'
            #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
            labs <- unique(d$CHR)
          }
          
          xmax = ceiling(max(d$pos) * 1.03)
          xmin = floor(max(d$pos) * -0.03)
          
          ########draw plot#######################
          
          par(mar=mar)
          def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                           xlim=c(xmin,xmax), ylim=c(0,log_times*max(d$logp)),
                           xlab=xlabel,ylab="",mgp=c(HorLabDis,0,0),cex.lab=label_size)
          
          dotargs <- list(NULL)
          do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
          axis(1, at=ticks, labels=labs,lwd=CoorLwd,tck=TckLen,mgp=c(2.5,HorTckDis,0.5),cex.axis=TckLwd)
          
          
         
          suppressWarnings(axis(2, at=seq(0,log_times*max(d$logp),ceiling(log_times*max(d$logp)/5)),lwd=CoorLwd,tck=TckLen,mgp=c(2.2,VerTckDis,0),cex.axis=TckLwd))
          mtext(expression(-log[10]('P-value')),side=2,line=VerLabDis,cex=label_size,font=1)
          
          # Create a vector of alternatiting colors
          col=rep(col, max(d$CHR))
          # Add points to the plot
          if (nchr==1) {
            with(d, points(pos, logp, pch=20, col=col[1]))
          } else {
            # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
            icol=1
            for (i in unique(d$index)) {
              with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20))
              icol=icol+1
            }
          }
          d.highlight=d[which(d$SNP %in% highlight), ]
          highlight_LOD<-as.numeric(data_fin_mark_order[,3])
          d.highlight<-as.data.frame(cbind(d.highlight,highlight_LOD))
          
          ################################
          par(new=T)
          def_args <- list(xaxt='n', yaxt='n',bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                           xlim=c(xmin,xmax), ylim=c(0,LOD_times*max(highlight_LOD)),xlab="",ylab="")
          dotargs <- list(NULL)
          do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
          suppressWarnings(axis(4,mgp=c(1.4,VerTckDis,0),at=seq(0,LOD_times*max(highlight_LOD),ceiling(LOD_times*max(highlight_LOD)/5)),col="magenta",col.ticks="magenta",col.axis="magenta",lwd=CoorLwd,tck=TckLen,cex.axis=TckLwd))
          mtext("LOD score",side=4,line=VerLabDis,cex=label_size,font=1,col="magenta")
          abline(h=lodline,col="gray25",lty=2,lwd=2)
          peach_colors<-c("magenta","deepskyblue2")
          col_pos<-list(NULL)
          method_num<-sort(unique(data_fin_mark_order[,4]))
          
          if(max(unique(data_fin[,3]))<max(unique(data_fin_mark_order[,4]))){
            col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
            col_pos[[2]]<-which(data_fin_mark_order[,4]!=max(method_num))
          }else{
            if(length(unique(data_fin[,3]))==1){
              col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
            }else{
              col_pos[[1]]<-1:nrow(data_fin_mark_order)
              
            }
          }
          if(length(col_pos)>1&&length(col_pos[[2]])!=0){
            with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20))
            with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20,type="h",lty=2))
            with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20))
            with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,type="h",lty=2))
          }else{
            with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20))
            with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,type="h",lty=2))
          }
      
        }
        
        
        
        QQ_mrMLM<-function(data_in,mar=c(2.5,2.5,1,1),label_size=0.7,TckLen=-0.02,
                           CoorLwd=3,TckLwd=0.6,HorLabDis=1,HorTckDis=0.02,VerLabDis=1.1,
                           VerTckDis=0.3,P_stand=0.9){
          
          method<-unique(data_in[,3])
          data_method<-list(NULL)
          for(i in 1:length(method)){
            data_method[[i]]<-data_in[which(data_in[,3]==method[i]),]
          }
          logp_4method<-numeric()
          for(i in 1:length(method)){
            method_p<-data_method[[i]][,8]
            logp_4method<-cbind(logp_4method,method_p) 
          }
          logp_4method<-apply(logp_4method,2,as.numeric)
          p_4method<-10^-logp_4method
          p_median<-apply(p_4method,1,median)
          locsub<-which(p_median==0)
          pmin<-min(p_median[p_median!=0])
          subvalue<-10^(1.1*log10(pmin))
          p_median[locsub]<-subvalue
          data_p<-as.matrix(p_median)
          p_value<-data_p
          pvalue<-matrix(p_value,,1)
          observed<-sort(pvalue[,1])
          observed<-observed/2
          observed<-observed[which(observed!=0)]
          newobserved<-observed[which(observed<(0.9/2))]
          lobs<--(log10(newobserved))
          expected<-c(1:length(newobserved))
          lexp<--(log10(expected/(length(pvalue)+1)))
          par(mar=mar)
          suppressWarnings(plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P-value)'),
                                yaxt="n",ylab="",col="blue",pch=20,cex.lab=label_size,tck=TckLen,bty="l",lwd=CoorLwd,
                                lwd.ticks=CoorLwd,cex.axis=TckLwd,mgp=c(HorLabDis,HorTckDis,0)))
          suppressWarnings(axis(2, at=seq(0,max(lobs)),lwd=CoorLwd,tck=TckLen,mgp=c(2.2,VerTckDis,0),cex.axis=TckLwd))
          mtext(expression('Observed -log'[10]*'(P-value)'),side=2,line=VerLabDis,cex=label_size,font=1)
          abline(0,1,col="red")
          box(bty="l",lwd=CoorLwd)
        }
      }
      
      screen<-function(reMR,rawgen,gen_num,phe_num,ps_num){
        if(nrow(reMR)>=200){
          reMR4<-as.matrix(reMR[,4])
          datashuz1<-rawgen[-1,1]
          calculate_gene<-t(gen_num[which(datashuz1%in%reMR4),-c(1,2)])
          gene_shuzhi<-apply(calculate_gene,2,as.numeric)
          larsres<-lars(gene_shuzhi,phe_num,type = "lar",trace = FALSE,use.Gram=FALSE,max.steps=200)  
          X<-gene_shuzhi[,which(larsres$beta[nrow(larsres$beta),]!=0)]
          MR200<-reMR[which(larsres$beta[nrow(larsres$beta),]!=0),]
          z<-cbind(matrix(1,nrow(gene_shuzhi),1),ps_num) 
          u1<-try({sblgwas(z,phe_num,X,t = -4,max.iter = 200,min.err = 1e-8)},silent=TRUE)
          if('try-error' %in% class(u1)){
            u1<-try({sblgwas(z,phe_num,X,t = -2,max.iter = 200,min.err = 1e-8)},silent=TRUE)
          }
          reMRshai<-MR200[which(u1$blup$p_wald<=0.01),]
          ind1<-which(larsres$beta[nrow(larsres$beta),]!=0)
          indz<-ind1[which(u1$blup$p_wald<=0.01)]
        }else if(nrow(reMR)<200){
          reMR4<-as.matrix(reMR[,4])
          datashuz1<-rawgen[-1,1]
          calculate_gene<-t(gen_num[which(datashuz1%in%reMR4),-c(1,2)])
          gene_shuzhi<-apply(calculate_gene,2,as.numeric)
          X<-gene_shuzhi
          z<-cbind(matrix(1,nrow(gene_shuzhi),1),ps_num) 
          u1<-try({sblgwas(z,phe_num,X,t = -4,max.iter = 200,min.err = 1e-8)},silent=TRUE)
          if('try-error' %in% class(u1)){
            u1<-try({sblgwas(z,phe_num,X,t = -2,max.iter = 200,min.err = 1e-8)},silent=TRUE)
          }
          reMRshai<-reMR[which(u1$blup$p_wald<=0.01),]
          
          indz<-which(u1$blup$p_wald<=0.01)
        }
        reMR<-cbind(reMRshai[,1:12],reMR[1:nrow(reMRshai),13:14])
        result<-list(reMR,indz)
        return(result)
      }
      
      readraw<-ReadData()
      if(input$radiomrformat=="mrMLM numeric format"){
        Genformat<-1  
      }else if(input$radiomrformat=="mrMLM character format"){
        Genformat<-2
      }else if(input$radiomrformat=="Hapmap (TASSEL) format"){
        Genformat<-3
      }
      Plotformat1<-Plotformat;Plotformat2<-Plotformat
      PheName<-readraw$phename
      
      if(length(grep(":",input$trait))!=0){
        scope<-as.numeric(unlist(strsplit(input$trait,":")))
        trait<-as.numeric(scope[1]:scope[2])
      }else if(length(grep(",",input$trait))!=0){
        scope<-as.numeric(unlist(strsplit(input$trait,",")))
        trait<-c(scope)
      }else if(length(grep(":",input$trait))==0&&length(grep(",",input$trait))==0){
        trait<-as.numeric(input$trait)
      }
      print("Running in progress, please be patient...")
      
      withProgress(message = 'Running in progress', value = 0, {
        for (i in trait){
          InputData<-InputDataF(readraw,Genformat,input$Method,i,input$poptypeall)
          reMR<-NULL;reFMR<-NULL;reFME<-NULL;rePLA<-NULL;rePKW<-NULL;reISIS<-NULL
          re1MR<-NULL;re1FMR<-NULL;re1FME<-NULL;re1PLA<-NULL;re1PKW<-NULL;re1ISIS<-NULL
          remanMR<-NULL;reqqMR<-NULL;remanFMR<-NULL;reqqFMR<-NULL;remanFME<-NULL;reqqFME<-NULL;
          replPLA<-NULL;remanPKW<-NULL;reqqPKW<-NULL; replISIS<-NULL;metaresult<-NULL;result_output<-NULL
          
          method<-input$Method
          TRY1<-try({
            if("mrMLM"%in%method){
              outMR<-mrMLMFun(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.01,svrad,svmlod,Genformat,CLO=NULL)
              if(is.null(outMR$result2)==FALSE){
                me<-matrix("mrMLM",nrow(outMR$result2),1)
                tr<-matrix(i,nrow(outMR$result2),1)
                trna<-matrix(PheName[i,],nrow(outMR$result2),1)
                colnames(me)<-"Method"
                colnames(tr)<-"Trait ID"
                colnames(trna)<-"Trait name"
                reMR<-cbind(tr,trna,me,as.matrix(outMR$result2))
                if(nrow(reMR)>50){
                  reMR<-screen(reMR,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)[[1]]
                }
              }
              me1<-matrix("mrMLM",nrow(outMR$result1),1)
              tr1<-matrix(i,nrow(outMR$result1),1)
              tr1na<-matrix(PheName[i,],nrow(outMR$result1),1)
              colnames(me1)<-"Method"
              colnames(tr1)<-"Trait ID"
              colnames(tr1na)<-"Trait name"
              re1MR<-cbind(tr1,tr1na,me1,as.matrix(outMR$result1))
            }
          },silent=FALSE)  
          if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){
            TRY2<-try({
              if("FASTmrMLM"%in%method){
                outFMR<-FASTmrMLM(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.01,svrad,svmlod,Genformat,CLO=NULL)
                if(is.null(outFMR$result2)==FALSE){
                  me<-matrix("FASTmrMLM",nrow(outFMR$result2),1)
                  tr<-matrix(i,nrow(outFMR$result2),1)
                  trna<-matrix(PheName[i,],nrow(outFMR$result2),1)
                  colnames(me)<-"Method"
                  colnames(tr)<-"Trait ID"
                  colnames(trna)<-"Trait name"
                  reFMR<-cbind(tr,trna,me,as.matrix(outFMR$result2))
                  if(nrow(reFMR)>50){
                    reFMR<-screen(reFMR,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)[[1]]
                  }
                }
                me1<-matrix("FASTmrMLM",nrow(outFMR$result1),1)
                tr1<-matrix(i,nrow(outFMR$result1),1)
                tr1na<-matrix(PheName[i,],nrow(outFMR$result1),1)
                colnames(me1)<-"Method"
                colnames(tr1)<-"Trait ID"
                colnames(tr1na)<-"Trait name"
                re1FMR<-cbind(tr1,tr1na,me1,as.matrix(outFMR$result1))
              }
              
            },silent=FALSE)
          }
          if ('try-error' %in% class(TRY2)|| !('try-error' %in% class(TRY2))){
            TRY3<-try({
              if("FASTmrEMMA"%in%method){
                outFME<-FASTmrEMMA(InputData$doFME$gen,InputData$doFME$phe,InputData$doFME$outATCG,InputData$doFME$genRaw,InputData$doFME$kk,InputData$doFME$psmatrix,0.005,svmlod,Genformat,Likelihood,CLO=NULL)
                if(is.null(outFME$result2)==FALSE){
                  me<-matrix("FASTmrEMMA",nrow(outFME$result2),1)
                  tr<-matrix(i,nrow(outFME$result2),1)
                  trna<-matrix(PheName[i,],nrow(outFME$result2),1)
                  colnames(me)<-"Method"
                  colnames(tr)<-"Trait ID"
                  colnames(trna)<-"Trait name"
                  reFME<-cbind(tr,trna,me,as.matrix(outFME$result2))
                  if(nrow(reFME)>50){
                    reFME<-screen(reFME,InputData$doFME$genRaw,InputData$doFME$gen,InputData$doFME$phe,InputData$doFME$psmatrix)[[1]]
                  }
                }
                me1<-matrix("FASTmrEMMA",nrow(outFME$result1),1)
                tr1<-matrix(i,nrow(outFME$result1),1)
                tr1na<-matrix(PheName[i,],nrow(outFME$result1),1)
                colnames(me1)<-"Method"
                colnames(tr1)<-"Trait ID"
                colnames(tr1na)<-"Trait name"
                re1FME<-cbind(tr1,tr1na,me1,as.matrix(outFME$result1))
              }
            },silent=FALSE)
            
          }
          if ('try-error' %in% class(TRY3)|| !('try-error' %in% class(TRY3))){
            TRY4<-try({
              if("pLARmEB"%in%method){
                outPLA<-pLARmEB(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,svmlod,lars1,Genformat,Bootstrap,CLO = NULL)
                if(is.null(outPLA$result)==FALSE){
                  me<-matrix("pLARmEB",nrow(outPLA$result),1)
                  tr<-matrix(i,nrow(outPLA$result),1)
                  trna<-matrix(PheName[i,],nrow(outPLA$result),1)
                  colnames(me)<-"Method"
                  colnames(tr)<-"Trait ID"
                  colnames(trna)<-"Trait name"
                  rePLA<-cbind(tr,trna,me,as.matrix(outPLA$result))
                  replPLA<-outPLA$plot
                  if(nrow(rePLA)>50){
                    rePLAQ<-screen(rePLA,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)
                    rePLA<-rePLAQ[[1]]
                  }
                }
              }
            },silent=FALSE)
            
          }
          if ('try-error' %in% class(TRY4)|| !('try-error' %in% class(TRY4))){
            TRY5<-try({
              if("pKWmEB"%in%method){
                outPKW<-pKWmEB(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.05,svmlod,Genformat,CLO=NULL)
                if(is.null(outPKW$result2)==FALSE){
                  me<-matrix("pKWmEB",nrow(outPKW$result2),1)
                  tr<-matrix(i,nrow(outPKW$result2),1)
                  trna<-matrix(PheName[i,],nrow(outPKW$result2),1)
                  colnames(me)<-"Method"
                  colnames(tr)<-"Trait ID"
                  colnames(trna)<-"Trait name"
                  rePKW<-cbind(tr,trna,me,as.matrix(outPKW$result2))
                  if(nrow(rePKW)>50){
                    rePKW<-screen(rePKW,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)[[1]]
                  }
                }
                me1<-matrix("pKWmEB",nrow(outPKW$result1),1)
                tr1<-matrix(i,nrow(outPKW$result1),1)
                tr1na<-matrix(PheName[i,],nrow(outPKW$result1),1)
                colnames(me1)<-"Method"
                colnames(tr1)<-"Trait ID"
                colnames(tr1na)<-"Trait name"
                re1PKW<-cbind(tr1,tr1na,me1,as.matrix(outPKW$result1))
              }
            },silent=FALSE)
          }
          if ('try-error' %in% class(TRY5)|| !('try-error' %in% class(TRY5))){
            TRY6<-try({
              if("ISIS EM-BLASSO"%in%method){
                outISIS<-ISIS(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,0.01,svmlod,Genformat,CLO=NULL)
                if(is.null(outISIS$result)==FALSE){
                  me<-matrix("ISIS EM-BLASSO",nrow(outISIS$result),1)
                  tr<-matrix(i,nrow(outISIS$result),1)
                  trna<-matrix(PheName[i,],nrow(outISIS$result),1)
                  colnames(me)<-"Method"
                  colnames(tr)<-"Trait ID"
                  colnames(trna)<-"Trait name"
                  reISIS<-cbind(tr,trna,me,as.matrix(outISIS$result))
                  replISIS<-outISIS$plot
                  if(nrow(reISIS)>50){
                    reISISQ<-screen(reISIS,InputData$doMR$genRaw,InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$psmatrix)
                    reISIS<-reISISQ[[1]]
                  }
                }
              }
            },silent=FALSE)
          }
          if ('try-error' %in% class(TRY6)|| !('try-error' %in% class(TRY6))){
            TRY7<-try({
              output1qq<-list(re1MR,re1FMR,re1FME,re1PKW)
              output1q<-do.call(rbind,output1qq)
              
            if(isFALSE(all(lengths(output1qq)==0))){
              eff<-numeric()
              logp<-numeric()
              for(bb in c(which(lengths(output1qq)!=0))){
                eff_every<-as.matrix(output1qq[[bb]][,7])
                colnames(eff_every)<-colnames(output1qq[[bb]])[7]
                eff<-cbind(eff,eff_every)
                  
                logp_every<-as.matrix(output1qq[[bb]][,8])
                colnames(logp_every)<-colnames(output1qq[[bb]])[8]
                logp<-cbind(logp,logp_every)
              }
              gencode1<-as.matrix(output1qq[[which(lengths(output1qq)!=0)[1]]][,9])
              colnames(gencode1)<-colnames(output1q)[[9]]
              
              output1<-cbind(output1qq[[which(lengths(output1qq)!=0)[1]]][,c(1,2,4,5,6)],eff,logp,gencode1)
              if("SNP effect (pKWmEB)"%in%colnames(output1)){
                output1<-output1[,-c(which(colnames(output1)%in%"SNP effect (pKWmEB)"))] 
              }
            }else{
              output1<-output1q
            }
              write.table(output1,paste(dir,"/",i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
               
            },silent=FALSE)
          }
          
          if ('try-error' %in% class(TRY7)|| !('try-error' %in% class(TRY7))){
            TRY8<-try({
              
              output<-list(reMR,reFMR,reFME,rePLA,rePKW,reISIS)
              output<-do.call(rbind,output)
              write.table(output,paste(dir,"/",i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
               
            },silent=FALSE)
          }
          
          
          
          if ('try-error' %in% class(TRY8)|| !('try-error' %in% class(TRY8))){
            TRY9<-try({
              
              if(DrawPlot==TRUE){
                
                
              if(isFALSE(all(lengths(output1qq)==0))){
                manwidth<-28000;manhei<-7000;manwordre<-60;manfigurere<-600 
                qqwidth<-10000;qqhei<-10000;qqwordre<-60;qqfigurere<-600 
                
                if(Plotformat1=="*.png"){
                  png(paste(dir,"/",i,"_Manhattan plot.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=svmlod)
                  dev.off()
                  
                  png(paste(dir,"/",i,"_qq plot.png",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
                  QQ_mrMLM(data_in=as.matrix(output1q))
                  dev.off()
                  
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(dir,"/",i,"_Manhattan plot.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=svmlod)
                  dev.off()
                  
                  tiff(paste(dir,"/",i,"_qq plot.tiff",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
                  QQ_mrMLM(data_in=as.matrix(output1q))
                  dev.off()
                  
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(dir,"/",i,"_Manhattan plot.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),lodline=svmlod)
                  dev.off()
                  
                  jpeg(paste(dir,"/",i,"_qq plot.jpeg",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
                  QQ_mrMLM(data_in=as.matrix(output1q))
                  dev.off()
                  
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(dir,"/",i,"_Manhattan plot.pdf",sep=""),width=16,height=4,pointsize = 20)
                  manhattan_mrMLM(data_in=as.matrix(output1q),data_fin=as.matrix(output),CoorLwd=2,lodline=svmlod)
                  dev.off()
                  
                  pdf(paste(dir,"/",i,"_qq plot.pdf",sep=""),pointsize = 25)
                  QQ_mrMLM(data_in=as.matrix(output1q),CoorLwd=2)
                  dev.off()
                }
              
              }else{
                showModal(modalDialog(title = "Warning", strong("Draw plot need intermediate result of mrMLM, FASTmrMLM, FASTmrEMMA or pKWmEB!"), easyClose = TRUE))
              }   
            } 
              
            },silent=FALSE)
          }
          
           
          incProgress(1/(max(trait)-min(trait)+1), detail = paste("Doing part", i))
          Sys.sleep(0.1)
        }
      })
      
    }) 
    
    output$result<-renderUI(MR())
    
    ####################Draw plot###############################################
    
    manIn<-reactive({
      req(input$filemanIn)
      data_in<-as.matrix(fread(input$filemanIn$datapath))
    })
    
    manFin<-reactive({
      req(input$filemanFin)
      data_fin<-as.matrix(fread(input$filemanFin$datapath))
    })
    
    
    qqIn<-reactive({
      req(input$fileqqmr)
      data_in<-as.matrix(fread(input$fileqqmr$datapath))
    })
    
    
    mandraw<-eventReactive(input$drawman,{
      
      
      id <- showNotification("Calculation in progress, please be patient...", duration = NULL,type="message")
      
      manhattan_mrMLM2<-function(data_in,data_fin,mar=c(2.9,2.8,0.7,2.8),VerLabDis=1.5,
                                VerTckDis=0.4,label_size=0.8,CoorLwd=5,
                                TckLen=-0.03,TckLwd=0.7,log_times=2,LOD_times=1.2){
        
        ###########Data process#################
        ###########intermediate result
        cona<-colnames(data_in)
        logp_4method<-as.matrix(data_in[,which(substr(cona,1,10)=="'-log10(P)")])
        logp_4method<-apply(logp_4method,2,as.numeric)
        p_4method<-10^-logp_4method
        p_median<-apply(p_4method,1,median)
        locsub<-which(p_median==0)
        pmin<-min(p_median[p_median!=0])
        subvalue<-10^(1.1*log10(pmin))
        p_median[locsub]<-subvalue
        data_p<-as.matrix(p_median)
        data_num<-as.matrix(seq(1:length(p_median)))
        data_chr<-as.matrix(data_in[,4])
        data_pos<-as.matrix(data_in[,5])
        manresult<-cbind(data_chr,data_pos,data_p,data_num)
        manresult<-apply(manresult,2,as.numeric)
        colnames(manresult)<-c("Chromosome","BPnumber","P-value","SNPname")
        manresult<-as.data.frame(manresult)
        #######final result##################
        data_fin_method<-unique(data_fin[,3])
        data_fin_method_length<-1:length(unique(data_fin[,3]))
        for(r in 1:length(unique(data_fin[,3]))){
          data_fin[which(data_fin[,3]==data_fin_method[r]),3]<-r
        }
        data_fin_mark<-matrix(data_fin[,c(5,6,8,3)],,4)
        data_fin_mark<-matrix(apply(data_fin_mark,2,as.numeric),,4)
        data_fin_mark_chr<-matrix(data_fin_mark[order(data_fin_mark[,1]),],,4)
        data_fin_mark_order<-numeric()
        for(i in c(unique(data_fin_mark_chr[,1]))){
          data_fin_mark_erery_chr<-matrix(data_fin_mark_chr[which(data_fin_mark_chr[,1]==i),],,4)
          data_fin_mark_pos<-matrix(data_fin_mark_erery_chr[order(data_fin_mark_erery_chr[,2]),],,4)
          all_pos<-unique(data_fin_mark_pos[,2])
          all_pos_maxlod<-numeric()
          for(ii in 1:length(all_pos)){
            all_pos_every<-matrix(data_fin_mark_pos[which(data_fin_mark_pos[,2]==all_pos[ii]),],,4)
            lod_me<-median(all_pos_every[,3])
            all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,all_pos_every[1,4])
            if(nrow(all_pos_every)>=2){
              all_pos_every_median<-c(all_pos_every[1,1:2],lod_me,max(data_fin_mark[,4])+1)
            }
            all_pos_maxlod<-rbind(all_pos_maxlod,all_pos_every_median)
          }
          data_fin_mark_order<-rbind(data_fin_mark_order,all_pos_maxlod)
        }
        snpOfInterest<-numeric()
        for(i in c(unique(data_fin_mark_order[,1]))){
          manresult_chr<-manresult[which(manresult[,1]==i),]
          data_fin_mark_order_chr<-matrix(data_fin_mark_order[which(data_fin_mark_order[,1]==i),],,4)
          mark_loc<-manresult_chr[which(manresult_chr[,2]%in%data_fin_mark_order_chr[,2]),4]
          snpOfInterest<-c(snpOfInterest,mark_loc) 
        }
        bpnumber <- numeric()
        chrnum <- unique(manresult[,1])
        for(i in 1:length(chrnum))
        {
          bpnumber <- rbind(bpnumber,as.matrix(c(1:length(which(manresult[,1]==chrnum[i])))))
        }
        manresult2<-cbind(manresult[,1],bpnumber,manresult[,3:4])
        colnames(manresult2)<-c("Chromosome","BPnumber","P-value","SNPname")
        ##########prepare for data#############################
        x<-manresult2;col=c("lightgreen","lightskyblue");logp=TRUE
        chr = "Chromosome";bp ="BPnumber";p ="P-value";snp="SNPname";
        highlight<-snpOfInterest
        CHR=BP=P=index=NULL
        d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
        if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
        d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
        d <- d[order(d$CHR, d$BP), ]
        if (logp) {
          d$logp <- -log10(d$P)
        } else {
          d$logp <- d$P
        }
        d$pos=NA
        d$index=NA
        ind = 0
        for (i in unique(d$CHR)){
          ind = ind + 1
          d[d$CHR==i,]$index = ind
        }
        
        nchr = length(unique(d$CHR))
        if (nchr==1) { ## For a single chromosome
          ## Uncomment the next two linex to plot single chr results in Mb
          #options(scipen=999)
          #d$pos=d$BP/1e6
          d$pos=d$BP
          ticks=floor(length(d$pos))/2+1
          xlabel = paste('Chromosome',unique(d$CHR),'position')
          labs = ticks
        } else { ## For multiple chromosomes
          lastbase=0
          ticks=NULL
          for (i in unique(d$index)) {
            if (i==1) {
              d[d$index==i, ]$pos=d[d$index==i, ]$BP
            } else {
              lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
              d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
            }
            # Old way: assumes SNPs evenly distributed
            # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
            # New way: doesn't make that assumption
            ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
          }
          xlabel = 'Chromosomes'
          #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
          labs <- unique(d$CHR)
        }
        
        xmax = ceiling(max(d$pos) * 1.03)
        xmin = floor(max(d$pos) * -0.03)
        
        ########draw plot#######################
        
        par(mar=mar)
        def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                         xlim=c(xmin,xmax), ylim=c(0,log_times*max(d$logp)),
                         xlab=xlabel,ylab="",mgp=c(VerLabDis,0,0),cex.lab=label_size)
        
        dotargs <- list(NULL)
        do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
        axis(1, at=ticks, labels=labs,lwd=CoorLwd,tck=TckLen,mgp=c(VerLabDis+1,VerTckDis-0.2,0.5),cex.axis=TckLwd)
        
        suppressWarnings(axis(2, at=seq(0,log_times*max(d$logp),ceiling(log_times*max(d$logp)/5)),lwd=CoorLwd,tck=TckLen,mgp=c(2.2,VerTckDis,0),cex.axis=TckLwd))
        mtext(expression(-log[10]('P-value')),side=2,line=VerLabDis,cex=label_size,font=1)
        
        # Create a vector of alternatiting colors
        col=rep(col, max(d$CHR))
        # Add points to the plot
        if (nchr==1) {
          with(d, points(pos, logp, pch=20, col=col[1]))
        } else {
          # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
          icol=1
          for (i in unique(d$index)) {
            with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20))
            icol=icol+1
          }
        }
        d.highlight=d[which(d$SNP %in% highlight), ]
        highlight_LOD<-as.numeric(data_fin_mark_order[,3])
        d.highlight<-as.data.frame(cbind(d.highlight,highlight_LOD))
        
        ################################
        par(new=T)
        
        def_args <- list(xaxt='n', yaxt='n',bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                         xlim=c(xmin,xmax), ylim=c(0,LOD_times*max(highlight_LOD)),xlab="",ylab="")
        dotargs <- list(NULL)
        do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
        suppressWarnings(axis(4,mgp=c(1.4,VerTckDis,0),at=seq(0,LOD_times*max(highlight_LOD),ceiling(LOD_times*max(highlight_LOD)/5)),col="magenta",col.ticks="magenta",col.axis="magenta",lwd=CoorLwd,tck=TckLen,cex.axis=TckLwd))
        mtext("LOD score",side=4,line=VerLabDis,cex=label_size,font=1,col="magenta")
        abline(h=3,col="gray25",lty=2,lwd=2)
        peach_colors<-c("magenta","deepskyblue2")
        col_pos<-list(NULL)
        method_num<-sort(unique(data_fin_mark_order[,4]))
        
        if(max(unique(data_fin[,3]))<max(unique(data_fin_mark_order[,4]))){
          col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
          col_pos[[2]]<-which(data_fin_mark_order[,4]!=max(method_num))
        }else{
          if(length(unique(data_fin[,3]))==1){
            col_pos[[1]]<-which(data_fin_mark_order[,4]==max(method_num))
          }else{
            col_pos[[1]]<-1:nrow(data_fin_mark_order)
            
          }
        }
        if(length(col_pos)>1&&length(col_pos[[2]])!=0){
          with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20))
          with(d.highlight, points(pos[col_pos[[2]]], highlight_LOD[col_pos[[2]]], col=peach_colors[2], pch=20,type="h",lty=2))
          with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20))
          with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,type="h",lty=2))
        }else{
          with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20))
          with(d.highlight, points(pos[col_pos[[1]]], highlight_LOD[col_pos[[1]]], col=peach_colors[1], pch=20,type="h",lty=2))
        }
        
      }
      
      dirplot<-input$Saveplot
      
      manwidth<-input$manwidthmr;manhei<-input$manheimr;
      manwordre<-input$manwordremr;manfigurere<-input$manfigureremr
      
      VerLabDis<-as.numeric(input$VerLabDis);VerTckDis<-as.numeric(input$VerTckDis)
      label_size<-as.numeric(input$labelsize);CoorLwd<-as.numeric(input$CoorLwd)
      TckLen<-as.numeric(input$TckLen);TckLwd<-as.numeric(input$TckLwd)
      log_times<-as.numeric(input$logtimes);LOD_times<-as.numeric(input$LODtimes)
      
      if(input$markgene==TRUE){
        
        if(length(grep(",",input$posX))!=0&length(grep(",",input$posY))!=0&length(grep(",",input$Genename))!=0){
        pos_x<-as.numeric(unlist(strsplit(input$posX,",")))
        pos_y<-as.numeric(unlist(strsplit(input$posY,",")))
        Genename<-as.character(unlist(strsplit(input$Genename,","))) 
        pos_xx<-c(pos_x)
        pos_yy<-c(pos_y)
        GeneName<-c(Genename)
        }else{
        pos_xx<-as.numeric(input$posX)
        pos_yy<-as.numeric(input$posY)
        GeneName<-as.character(input$Genename)  
        }
      }
      
      if(input$plmanformat=="*.png"){
        png(paste(dirplot,"/","Manhattan plot.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
        manhattan_mrMLM2(data_in=manIn(),data_fin=manFin(),mar=c(2.9,2.8,0.7,2.8),
                        VerLabDis,VerTckDis,label_size,CoorLwd,
                        TckLen,TckLwd,log_times,LOD_times)
        if(input$markgene==TRUE){
          text(c(pos_xx),c(pos_yy),c(GeneName),font=3,cex=label_size-0.3,col=input$colorname)
         }
        dev.off()
      }else if(input$plmanformat=="*.tiff"){
        tiff(paste(dirplot,"/","Manhattan plot.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
        
        manhattan_mrMLM2(data_in=manIn(),data_fin=manFin(),mar=c(2.9,2.8,0.7,2.8),
                        VerLabDis,VerTckDis,label_size,CoorLwd,
                        TckLen,TckLwd,log_times,LOD_times)
        if(input$markgene==TRUE){
          text(c(pos_xx),c(pos_yy),c(GeneName),font=3,cex=label_size-0.3,col=input$colorname)
        }
        dev.off()
        
      }else if(input$plmanformat=="*.jpeg"){
        jpeg(paste(dirplot,"/","Manhattan plot.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
        
        manhattan_mrMLM2(data_in=manIn(),data_fin=manFin(),mar=c(2.9,2.8,0.7,2.8),
                        VerLabDis,VerTckDis,label_size,CoorLwd,
                        TckLen,TckLwd,log_times,LOD_times)
        if(input$markgene==TRUE){
          text(c(pos_xx),c(pos_yy),c(GeneName),font=3,cex=label_size-0.3,col=input$colorname)
        }
        dev.off()
        
      }else if(input$plmanformat=="*.pdf"){
        pdf(paste(dirplot,"/","Manhattan plot.pdf",sep=""),width=as.numeric(input$manwidthpdf),height=as.numeric(input$manheipdf),pointsize=as.numeric(input$manpoipdf))
        
        manhattan_mrMLM2(data_in=manIn(),data_fin=manFin(),mar=c(2.9,2.8,0.7,2.8),
                        VerLabDis,VerTckDis,label_size,CoorLwd,
                        TckLen,TckLwd,log_times,LOD_times)
        if(input$markgene==TRUE){
          text(c(pos_xx),c(pos_yy),c(GeneName),font=3,cex=label_size-0.3,col=input$colorname)
        }
        dev.off()
      }
      
      removeNotification(id)   
      
    })
    
    
    
    qqdraw<-eventReactive(input$drawqq,{
      
      id2 <- showNotification("Calculation in progress, please be patient...", duration = NULL,type="message")
      
      dirqq<-input$Saveplotqq
      
      QQ_mrMLM<-function(data_in,mar=c(2.5,2.5,1,1),label_size=0.7,TckLen=-0.02,
                         CoorLwd=3,TckLwd=0.6,HorLabDis=1,HorTckDis=0.02,VerLabDis=1.1,
                         VerTckDis=0.3,P_stand=0.9){
        
        cona<-colnames(data_in)
        logp_4method<-as.matrix(data_in[,which(substr(cona,1,10)=="'-log10(P)")])
        logp_4method<-apply(logp_4method,2,as.numeric)
        p_4method<-10^-logp_4method
        p_median<-apply(p_4method,1,median)
        locsub<-which(p_median==0)
        pmin<-min(p_median[p_median!=0])
        subvalue<-10^(1.1*log10(pmin))
        p_median[locsub]<-subvalue
        data_p<-as.matrix(p_median)
        p_value<-data_p
        pvalue<-matrix(p_value,,1)
        observed<-sort(pvalue[,1])
        observed<-observed/2
        observed<-observed[which(observed!=0)]
        newobserved<-observed[which(observed<(0.9/2))]
        lobs<--(log10(newobserved))
        expected<-c(1:length(newobserved))
        lexp<--(log10(expected/(length(pvalue)+1)))
        par(mar=mar)
        suppressWarnings(plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P-value)'),
                              yaxt="n",ylab="",col="blue",pch=20,cex.lab=label_size,tck=TckLen,bty="l",lwd=CoorLwd,
                              lwd.ticks=CoorLwd,cex.axis=TckLwd,mgp=c(HorLabDis,HorTckDis,0)))
        suppressWarnings(axis(2, at=seq(0,max(lobs)),lwd=CoorLwd,tck=TckLen,mgp=c(2.2,VerTckDis,0),cex.axis=TckLwd))
        mtext(expression('Observed -log'[10]*'(P-value)'),side=2,line=VerLabDis,cex=label_size,font=1)
        abline(0,1,col="red")
        box(bty="l",lwd=CoorLwd)
      }
      qqwidth<-input$qqwidthmr;qqhei<-input$qqheimr;qqwordre<-input$qqwordremr;qqfigurere<-input$qqfigureremr
      
      VerLabDis<-as.numeric(input$VerLabDisqq);HorLabDis<-as.numeric(input$HorLabDisqq)
      HorTckDis<-as.numeric(input$HorTckDisqq);VerTckDis<-as.numeric(input$VerTckDisqq)
      label_size<-as.numeric(input$labelsizeqq);CoorLwd<-as.numeric(input$CoorLwdqq)
      TckLen<-as.numeric(input$TckLenqq);TckLwd<-as.numeric(input$TckLwdqq)
      P_stand<-as.numeric(input$pstand);
      
      if(input$plqqformat=="*.png"){
        png(paste(dirqq,"/","qq plot.png",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
        QQ_mrMLM(data_in=qqIn(),mar=c(2.5,2.5,1,1),label_size,TckLen,CoorLwd,TckLwd,HorLabDis,HorTckDis,VerLabDis,VerTckDis,P_stand)
        dev.off()
      }else if(input$plqqformat=="*.tiff"){
        tiff(paste(dirqq,"/","qq plot.tiff",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
        QQ_mrMLM(data_in=qqIn(),mar=c(2.5,2.5,1,1),label_size,TckLen,CoorLwd,TckLwd,HorLabDis,HorTckDis,VerLabDis,VerTckDis,P_stand)
        dev.off()
      }else if(input$plqqformat=="*.jpeg"){
        jpeg(paste(dirqq,"/","qq plot.jpeg",sep=""),width=as.numeric(qqwidth), height=as.numeric(qqhei), units= "px", pointsize =as.numeric(qqwordre),res=as.numeric(qqfigurere))
        QQ_mrMLM(data_in=qqIn(),mar=c(2.5,2.5,1,1),label_size,TckLen,CoorLwd,TckLwd,HorLabDis,HorTckDis,VerLabDis,VerTckDis,P_stand)
        dev.off()
      }else if(input$plqqformat=="*.pdf"){
        pdf(paste(dirqq,"/","qq plot.pdf",sep=""),width=as.numeric(input$qqwidthpdf),height=as.numeric(input$qqheipdf),pointsize=as.numeric(input$qqpoipdf))
        QQ_mrMLM(data_in=qqIn(),mar=c(2.5,2.5,1,1),label_size,TckLen,CoorLwd,TckLwd,HorLabDis,HorTckDis,VerLabDis,VerTckDis,P_stand)
        dev.off()
      }
      
      removeNotification(id2)
    })
    
    output$manout<-renderPlot(mandraw())
    output$qqout<-renderPlot(qqdraw())
 
  })
  
  ui<-(
    fluidPage(
      
      shinyjs::useShinyjs(),
      radioButtons("layout", "",
                   choiceNames = list(
                     strong("mrMLM"),
                     strong("Start")
                     
                   ),
                   choiceValues = list(
                     "mrMLM", "Start"
                   ),inline=TRUE),
      uiOutput("general_ui")
    )
  )
  
  
  
  shinyApp(ui<-ui,server<-server)
  
  
}