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
                           naturally occurring genetic variance among commercial inbred lines of maize (Zea mays L.).
                           Genetics 2005;169:2267-2275"),
                        h4("2. Wang SB, Feng JY, Ren WL, Huang B, Zhou L, Wen YJ, Zhang J, Jim M Dunwell, Xu S*, Zhang YM*.
                           Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear
                           model methodology. Scientific Reports 2016;6:19444. doi:10.1038/srep19444 (mrMLM)"),
                        h4("3. Tamba CL, Ni YL, Zhang YM*. Iterative sure independence screening EM-Bayesian LASSO algorithm for
                           multi-locus genome-wide association studies. PLoS Computational Biology 2017;13(1):e1005357. doi:10.1371/
                           journal.pcbi.1005357 (ISIS EM-BLASSO)"),
                        h4("4. Zhang J, Feng JY, Ni YL, Wen YJ, Niu Y, Tamba CL, Yue C, Song QJ, Zhang YM*. pLARmEB: 
                           integration of least angle regression with empirical Bayes for multi-locus genome-wide association
                           studies. Heredity 2017;118(6):517-524. doi:10.1038/hdy.2017.8 (pLARmEB)"),
                         h4("5. Ren WL, Wen YJ, Jim M Dunwell, Zhang YM*. pKWmEB: integration of Kruskal-Wallis test with
                           empirical Bayes under polygenic background control for multi-locus genome-wide association study. 
                           Heredity 2018;120(3):208-218. https://doi.org/10.1038/s41437-017-0007-4 (pKWmEB)"),
                         h4("6. Wen YJ, Zhang H, Ni YL, Huang B, Zhang J, Feng JY, Wang SB, Jim M Dunwell, Zhang YM*, Wu R*.
						                Methodological implementation of mixed linear models in multi-locus genome-wide association studies. 
                           Briefings in Bioinformatics 2018;19(4):700-712 doi:10.1093/bib/bbw145 (FASTmrEMMA)"),
                         h4("7. Tamba CL, Zhang YM. A fast mrMLM algorithm for multi-locus genome-wide association studies.
                           bioRxiv;preprint first posted online Jun. 7, 2018;doi: https://doi.org/10.1101/341784. (FASTmrMLM)"),
                         h4("8. Zhang Ya-Wen, Tamba Cox Lwaka, Wen Yang-Jun, Li Pei, Ren Wen-Long, Ni Yuan-Li, Gao Jun, Zhang Yuan-Ming*. 
                           mrMLM v4.0: An R platform for multi-locus genome-wide association studies. Genomics, Proteomics & Bioinformatics
						               2019;resubmission"),
                         br(),
                         br(),
                        
                        h4("Authors: Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming"),
                        
                        h4("Maintainer: Zhang Yuan-Ming (soyzhang at mail.hzau.edu.cn)"), 
                        
                        h4("mrMLM.GUI version 4.0, Realeased October 2019"),
                        
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
                                               radioButtons("Resolution", "Plot resolution (All methods)", choices = c("High", "Low")),
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
                       column(12,
                              radioButtons("Precisionman", "Resolution of plot", c("General resolution", "High resolution", "Set by yourself"),inline = TRUE)
                       ),
                       column(3,
                              textInput("manwidthmr", label = "Figure width (px):",value = "960"),
                              selectInput(inputId = "manchcolour1mr",label = "Chromosome colour (odd):",
                                          choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "blue"),
                              fileInput("filemanmr", "File to draw Manhattan plot",multiple = TRUE,accept = ".csv")
                       ),
                       column(3,
                              textInput("manheimr", label = "Figure height (px)",value = "600"),
                              selectInput(inputId = "manchcolour2mr",label = "Chromosome colour (even):",
                                          choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "red")
                              
                       ),
                       column(3,
                              textInput("manwordremr", label = "Word resolution (1/72 inch, ppi):",value = "20"),
                              textInput("mancrimr", label = "Critical value for Manhattan Plot:",value = "3")
                       ),
                       
                       column(3,
                              textInput("manfigureremr", label = "Figure resolution (ppi):",value = "72"),
                              selectInput(inputId = "plotmethodman",label = "Method selection",
                                          choices = c("mrMLM","FASTmrMLM","FASTmrEMMA","pKWmEB"),selected = "mrMLM")
                       ),
                       column(12,
                              radioButtons("plmanformat", "Manhattan plot format", c("*.png", "*.tiff", "*.jpeg","*.pdf"),inline = TRUE)
                              
                       ),
                       
                       column(12,
                              downloadButton("downloadmanplot", "Save Manhattan plot")
                       ),
                       
                       column(12,
                              plotOutput("mplmr")
                              
                       )
                     )
            ),
            
            tabPanel("QQ Plot",id="qqp",
                     
                     fluidRow(
                       
                       column(12,
                              radioButtons("Precisionqq", "Resolution of plot", c("General resolution", "High resolution", "Set by yourself"),inline = TRUE) 
                       ),
                       
                       column(3,
                              textInput("qqwidthmr", label = "Figure width (px):",value = "960"),
                              selectInput(inputId = "qqchcolour1mr",label = "Point colour:",
                                          choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "red"),
                              fileInput("fileqqmr", "File to draw QQ plot",multiple = TRUE,accept = ".csv")
                       ),
                       column(3,
                              textInput("qqheimr", label = "Figure height (px)",value = "600"),
                              selectInput(inputId = "qqchcolour2mr",label = "Line colour:",
                                          choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "black")
                              
                       ),
                       column(3,
                              textInput("qqwordremr", label = "Word resolution (1/72 inch, ppi):",value = "20"),
                              selectInput(inputId = "plotmethodqq",label = "Method selection",
                                          choices = c("mrMLM","FASTmrMLM","FASTmrEMMA","pKWmEB"),selected = "mrMLM")
                       ),
                       
                       column(3,
                              textInput("qqfigureremr", label = "Figure resolution (ppi):",value = "72")
                       ),
                       column(12,
                              radioButtons("plqqformat", "QQ plot format", c("*.png", "*.tiff", "*.jpeg","*.pdf"),inline = TRUE)
                       ),
                       
                       column(12,
                              downloadButton("downloadqqplot", "Save QQ plot")
                       ),
                       column(12,
                              plotOutput("qplmr")                                
                       )
                     ) 
            ),
            
            tabPanel("Plot of LOD Score against Genome position",id="lodp",
                     
                     fluidRow(
                       
                       column(12,
                              radioButtons("Precisionlod", "Resolution of plot", c("General resolution", "High resolution", "Set by yourself"),inline = TRUE) 
                       ),
                       
                       column(4,
                              textInput("ppwidthPLA", label = "Figure width (px):",value = "960"),
                              textInput("ppfigurerePLA", label = "Figure resolution (ppi):",value = "72"),
                              fileInput("fileppPLA", "File to draw plot",multiple = TRUE,accept = ".csv")
                       ),
                       column(4,
                              textInput("ppheiPLA", label = "Figure height (px)",value = "240"),
                              selectInput(inputId = "ppchcolour1PLA",label = "Lod line colour:",
                                          choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "red")
                              
                       ),
                       column(4,
                              textInput("ppwordrePLA", label = "Word resolution (1/72 inch, ppi):",value = "12"),
                              selectInput(inputId = "plotmethodlod",label = "Method selection",
                                          choices = c("pLARmEB","ISIS EM-BLASSO"),selected = "pLARmEB")
                       ),
                       column(12,
                              radioButtons("plppformatPLA", "Plot format", c("*.png", "*.tiff", "*.jpeg","*.pdf"),inline = TRUE)
                              
                       ),
                       
                       column(12,
                              downloadButton("downloadppplotPLA", "Save plot")
                       ),                         
                       
                       column(12,
                              plotOutput("qplPLA")
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
          kkPre<-as.matrix(kkRaw[-1,-1])
          nameKin<-as.matrix(kkRaw[-1,1])
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
    
    upda1<-observeEvent(input$Precisionman,{
      
      if(input$Precisionman=="High resolution"){
        widthman<-10000
        heightman<-6000
        wordman<-30
        figureman<-300
      }else if(input$Precisionman=="General resolution"){
        widthman<-960
        heightman<-600
        wordman<-20
        figureman<-72
      }else if(input$Precisionman=="Set by yourself"){
        widthman<-0
        heightman<-0
        wordman<-0
        figureman<-0
      }
      
      updateTextInput(session, "manwidthmr", value=widthman)
      updateTextInput(session, "manheimr", value=heightman)
      updateTextInput(session, "manwordremr", value=wordman)
      updateTextInput(session, "manfigureremr", value=figureman)
      
    }) 
    
    
    upda2<-observeEvent(input$Precisionqq,{   
      
      if(input$Precisionqq=="High resolution"){
        widthqq<-10000
        heightqq<-6000
        wordqq<-30
        figureqq<-300
      }else if(input$Precisionqq=="General resolution"){
        widthqq<-960
        heightqq<-600
        wordqq<-20
        figureqq<-72
        
      }else if(input$Precisionqq=="Set by yourself"){
        widthqq<-0
        heightqq<-0
        wordqq<-0
        figureqq<-0
      }
      updateTextInput(session, "qqwidthmr", value=widthqq)
      updateTextInput(session, "qqheimr", value=heightqq)
      updateTextInput(session, "qqwordremr", value=wordqq)
      updateTextInput(session, "qqfigureremr", value=figureqq)
    }) 
    
    
    upda3<-observeEvent(input$Precisionlod,{   
      if(input$Precisionlod=="High resolution"){
        widthlod<-10000
        heightlod<-6000
        wordlod<-30
        figurelod<-600
      }else if(input$Precisionlod=="General resolution"){
        widthlod<-960
        heightlod<-240
        wordlod<-12
        figurelod<-72
        
      }else if(input$Precisionlod=="Set by yourself"){
        widthlod<-0
        heightlod<-0
        wordlod<-0
        figurelod<-0
      }
      updateTextInput(session, "ppwidthPLA", value=widthlod)
      updateTextInput(session, "ppheiPLA", value=heightlod)
      updateTextInput(session, "ppwordrePLA", value=wordlod)
      updateTextInput(session, "ppfigurerePLA", value=figurelod)
    })
    
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
        Plot<-function(plotresult=NULL,color1=NULL,color2=NULL,p_stand=NULL,method=NULL,type=NULL){
          Manhattan<-function(plotresult,color1,color2){
            parms<-as.data.frame(plotresult)
            mannewp<-as.numeric(parms[1,5])
            svgwline<-round(-log10(mannewp),4)  
            standline<-svgwline
            manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P-value",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = standline)
          }
          
          QQplot1<-function(plotresult,p_stand,color1,color2){
            p_value<-as.matrix(plotresult)
            pvalue<-matrix(p_value,,1)
            observed<-sort(pvalue[,1])
            observed<-observed/2
            observed<-observed[which(observed!=0)]
            newobserved<-observed[which(observed<(p_stand/2))]
            lobs<--(log10(newobserved))
            expected<-c(1:length(newobserved))
            lexp<--(log10(expected/(length(pvalue)+1)))
            plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color2)
            abline(0,1,col=color1)
          }
          QQplot2<-function(plotresult,color1,color2){
            ress1<-as.data.frame(plotresult)
            pvalue<-as.matrix(ress1)
            ps<-pvalue[,1]
            obs.x<-sort(ps)
            newobs.x<-obs.x[obs.x<1]
            n<-length(newobs.x)
            es<-(1:n)/(n+1)
            x<--log10(es)
            y<--log10(newobs.x)
            y<-y-0.3
            plot(x,y,xlim=c(0.3,max(x)),ylim=c(0.3,max(y)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color2)
            abline(0,1,col=color1)
            
          }
          LOD<-function(fileplot=NULL,color1,method=NULL){
            data<-as.matrix(plotresult)
            data<-as.data.frame(data,stringsAsFactors = F)
            gen<-data[,1:2]
            resulty<-data[,3:5]
            resultkq<-as.matrix(resulty)
            resultk<-which(resultkq=="",arr.ind = TRUE)
            resultq<-resulty[1:(resultk[1]-1),] 
            
            if(nrow(resultq)>1){
              result<-resultq
            }else{
              result<-t(as.matrix(resultq))
            }
            
            galaxyy<-as.data.frame(result)
            galaxyy<-sapply(galaxyy,as.numeric)
            chr_pos <- gen[,1:2]
            chr_pos<-sapply(chr_pos,as.numeric)
            
            chr_num <- length(unique(chr_pos[,1]))
            chr <- matrix(0,chr_num,1)
            pos <- matrix(0,chr_num,1)
            for(i in 1:chr_num)
            {
              
              temp <- numeric()
              temp <- length(which(chr_pos[,1]==i))
              if(i==1)
              {
                pos[i] <- temp
                chr[i] <- chr_pos[pos[i],2]
              }else{
                pos[i] <- pos[i-1] + temp
                chr[i] <- chr_pos[pos[i],2]
              }
            }
            
            pos_acc <- matrix(0,chr_num,1)
            for(i in 1:chr_num)
            {
              if(i==1){
                pos_acc[i] <- chr[i]
              }else{
                pos_acc[i] <- pos_acc[i-1] + chr[i]
              }
            }
            
            newres_pos <- galaxyy[,2]
            res_sumpos <- pos_acc[galaxyy[which(galaxyy[,1]>1),1]-1] + galaxyy[which(galaxyy[,1]>1),2] 
            newres_pos[which(galaxyy[,1]>1)] <- res_sumpos 
            pospic<-c(newres_pos)
            lodpic<-c(galaxyy[,3])  
            mm<-round(max(pospic)/4000)
            mm<-as.numeric(format(mm,digits = 1,scientific = TRUE))
            pospicx<-pospic/mm
            if(pospicx[1]<20){
              pospicx[1]<-pospicx[1]+20
            }
            pos_acc1<-pos_acc/mm
            resdf1 <- data.frame(pospicx,lodpic)
            
            pp <- ggplot(data=resdf1, aes(x=pospicx, y=lodpic)) +
              geom_bar(stat="identity", width=0.5, fill="white", linetype="solid",color=color1)
            
            pp <- pp + geom_vline(xintercept=c(0,pos_acc1),linetype="dashed",alpha=0.2)
            pp <- pp  + scale_x_continuous(expand=c(0,0),limits=c(0,(pos_acc1[dim(pos_acc1)[1]]+100))) +
              scale_y_continuous(expand=c(0,0))
            pp <- pp + xlab(paste("Genome position (",mm,"bp)",sep = "")) + ylab("LOD score") + ggtitle("") + theme_classic()
            pp <- pp + theme(axis.title.y = element_text( vjust = 2,hjust=0.5,size = 14),
                             axis.title.x = element_text(vjust = -0.5,hjust=0.5,size = 14))
            
            pp <- pp + theme(panel.background = element_rect(fill = "white"))
            pp <- pp + theme(text=element_text(family="mono"))
            pp <- pp + theme(axis.line.y = element_line(colour = "black", linetype = "solid"),
                             axis.line.x = element_line(colour = "black", linetype = "solid"))
            print(pp)
            
          }  
          
          if(type=="Manhattan"){
            Manhattan(plotresult,color1,color2)
          }else if(type=="qq"){
            if(method=="FASTmrEMMA"){
              QQplot2(plotresult,color1,color2)  
            }else{
              QQplot1(plotresult,p_stand,color1,color2) 
            }
          }else if(type=="LOD"){
            LOD(plotresult,color1)
          }  
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
              remanMR<-outMR$Manhattan
              reqqMR<-outMR$QQ
              if(DrawPlot==TRUE){
                if(Resolution=="Low"){
                  manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
                }else if(Resolution=="High"){
                  manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300
                }
                if(Plotformat1=="*.png"){
                  png(paste(dir,"/",i,"_mrMLM_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(dir,"/",i,"_mrMLM_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(dir,"/",i,"_mrMLM_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(dir,"/",i,"_mrMLM_Manhattan.pdf",sep=""),width=10)
                }
                Plot(plotresult=remanMR,color1="red",color2="blue",0.95,method="mrMLM",type="Manhattan")
                dev.off()
                if(Plotformat2=="*.png"){
                  png(paste(dir,"/",i,"_mrMLM_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.tiff"){
                  tiff(paste(dir,"/",i,"_mrMLM_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.jpeg"){
                  jpeg(paste(dir,"/",i,"_mrMLM_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.pdf"){
                  pdf(paste(dir,"/",i,"_mrMLM_qq.pdf",sep=""),width=10)
                }
                Plot(plotresult=reqqMR,color1="red",color2="blue",0.95,method="mrMLM",type="qq")
                dev.off()
              }
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
                remanFMR<-outFMR$Manhattan
                reqqFMR<- outFMR$QQ
                if(DrawPlot==TRUE){
                  if(Resolution=="Low"){
                    manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
                  }else if(Resolution=="High"){
                    manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300
                  }
                  if(Plotformat1=="*.png"){
                    png(paste(dir,"/",i,"_FASTmrMLM_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.tiff"){
                    tiff(paste(dir,"/",i,"_FASTmrMLM_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_FASTmrMLM_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.pdf"){
                    pdf(paste(dir,"/",i,"_FASTmrMLM_Manhattan.pdf",sep=""),width=10)
                  }
                  Plot(plotresult=remanFMR,color1="red",color2="blue",0.95,method="FASTmrMLM",type="Manhattan")
                  dev.off()
                  if(Plotformat2=="*.png"){
                    png(paste(dir,"/",i,"_FASTmrMLM_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.tiff"){
                    tiff(paste(dir,"/",i,"_FASTmrMLM_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_FASTmrMLM_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.pdf"){
                    pdf(paste(dir,"/",i,"_FASTmrMLM_qq.pdf",sep=""),width=10)
                  }
                  Plot(plotresult=reqqFMR,color1="red",color2="blue",0.95,method="FASTmrMLM",type="qq")
                  dev.off()
                }
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
                remanFME<-outFME$Manhattan
                reqqFME<-outFME$QQ
                if(DrawPlot==TRUE){
                  if(Resolution=="Low"){
                    manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
                  }else if(Resolution=="High"){
                    manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300
                  }
                  if(Plotformat1=="*.png"){
                    png(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.tiff"){
                    tiff(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.pdf"){
                    pdf(paste(dir,"/",i,"_FASTmrEMMA_Manhattan.pdf",sep=""),width=10)
                  }
                  Plot(plotresult=remanFME,color1="red",color2="blue",0.95,method="FASTmrEMMA",type="Manhattan")
                  dev.off()
                  if(Plotformat2=="*.png"){
                    png(paste(dir,"/",i,"_FASTmrEMMA_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.tiff"){
                    tiff(paste(dir,"/",i,"_FASTmrEMMA_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_FASTmrEMMA_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.pdf"){
                    pdf(paste(dir,"/",i,"_FASTmrEMMA_qq.pdf",sep=""),width=10)
                  }
                  Plot(plotresult=reqqFME,color1="red",color2="blue",0.95,method="FASTmrEMMA",type="qq")
                  dev.off()
                }
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
                    plotid<-rePLAQ[[2]] 
                    loc_marker<-rbind(outPLA$plot[plotid,3:5],matrix("",nrow(outPLA$plot[,1:2])-length(plotid),3))
                    replPLA<-cbind(outPLA$plot[,1:2], loc_marker)
                  }
                }
                if(DrawPlot==TRUE){
                  if(Resolution=="Low"){
                    manwidth<-960;manhei<-240;manwordre<-12;manfigurere<-72
                  }else if(Resolution=="High"){
                    manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-600
                  }
                  if(Plotformat1=="*.png"){
                    png(paste(dir,"/",i,"_pLARmEB_LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.tiff"){
                    tiff(paste(dir,"/",i,"_pLARmEB_LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_pLARmEB_LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.pdf"){
                    pdf(paste(dir,"/",i,"_pLARmEB_LOD.pdf",sep=""),width=12)
                  }
                  Plot(plotresult=replPLA,color1="red",color2="blue",0.95,method="pLARmEB",type="LOD")
                  dev.off()
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
                remanPKW<-outPKW$Manhattan
                reqqPKW<-outPKW$QQ
                if(DrawPlot==TRUE){
                  if(Resolution=="Low"){
                    manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
                  }else if(Resolution=="High"){
                    manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-300
                  }
                  if(Plotformat1=="*.png"){
                    png(paste(dir,"/",i,"_pKWmEB_Manhattan.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.tiff"){
                    tiff(paste(dir,"/",i,"_pKWmEB_Manhattan.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_pKWmEB_Manhattan.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.pdf"){
                    pdf(paste(dir,"/",i,"_pKWmEB_Manhattan.pdf",sep=""),width=10)
                  }
                  Plot(plotresult=remanPKW,color1="red",color2="blue",0.95,method="pKWmEB",type="Manhattan")
                  dev.off()
                  if(Plotformat2=="*.png"){
                    png(paste(dir,"/",i,"_pKWmEB_qq.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.tiff"){
                    tiff(paste(dir,"/",i,"_pKWmEB_qq.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_pKWmEB_qq.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat2=="*.pdf"){
                    pdf(paste(dir,"/",i,"_pKWmEB_qq.pdf",sep=""),width=10)
                  }
                  Plot(plotresult=reqqPKW,color1="red",color2="blue",0.95,method="pKWmEB",type="qq")
                  dev.off()
                }
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
                    plotid<-reISISQ[[2]]  
                    loc_marker<-rbind(outISIS$plot[plotid,3:5],matrix("",nrow(outISIS$plot[,1:2])-length(plotid),3))
                    replISIS<-cbind(outISIS$plot[,1:2], loc_marker)
                  }
                }
                if(DrawPlot==TRUE){
                  if(Resolution=="Low"){
                    manwidth<-960;manhei<-240;manwordre<-12;manfigurere<-72
                  }else if(Resolution=="High"){
                    manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-600
                  }
                  if(Plotformat1=="*.png"){
                    png(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.tiff"){
                    tiff(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.jpeg"){
                    jpeg(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                  }else if(Plotformat1=="*.pdf"){
                    pdf(paste(dir,"/",i,"_ISIS EM-BLASSO_LOD.pdf",sep=""),width=12)
                  }
                  Plot(plotresult=replISIS,color1="red",color2="blue",0.95,method="ISIS EM-BLASSO",type="LOD")
                  dev.off()
                }
              }
            },silent=FALSE)
          }
          if ('try-error' %in% class(TRY6)|| !('try-error' %in% class(TRY6))){
            TRY7<-try({
              output1<-list(re1MR,re1FMR,re1FME,re1PKW)
              output1<-do.call(rbind,output1)
              output<-list(reMR,reFMR,reFME,rePLA,rePKW,reISIS)
              output<-do.call(rbind,output)
              write.table(output,paste(dir,"/",i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
              write.table(output1,paste(dir,"/",i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
              plotresult<-list(remanMR,reqqMR,remanFMR,reqqFMR,remanFME,reqqFME,replPLA,remanPKW,reqqPKW,replISIS)
              wb <- createWorkbook("Fred")
              addWorksheet(wb, "Manhattan mrMLM")
              addWorksheet(wb, "qq mrMLM")
              addWorksheet(wb, "Manhattan FASTmrMLM")
              addWorksheet(wb, "qq FASTmrMLM")
              addWorksheet(wb, "Manhattan FASTmrEMMA")
              addWorksheet(wb, "qq FASTmrEMMA")
              addWorksheet(wb, "Plot pLARmEB")
              addWorksheet(wb, "Manhattan pKWmEB")
              addWorksheet(wb, "qq pKWmEB")
              addWorksheet(wb, "Plot ISIS EM-BLASSO")
              writeData(wb, sheet = "Manhattan mrMLM", plotresult[[1]])
              writeData(wb, sheet = "qq mrMLM", plotresult[[2]])
              writeData(wb, sheet = "Manhattan FASTmrMLM", plotresult[[3]])
              writeData(wb, sheet = "qq FASTmrMLM", plotresult[[4]])
              writeData(wb, sheet = "Manhattan FASTmrEMMA", plotresult[[5]])
              writeData(wb, sheet = "qq FASTmrEMMA", plotresult[[6]])
              writeData(wb, sheet = "Plot pLARmEB", plotresult[[7]])
              writeData(wb, sheet = "Manhattan pKWmEB", plotresult[[8]])
              writeData(wb, sheet = "qq pKWmEB", plotresult[[9]])
              writeData(wb, sheet = "Plot ISIS EM-BLASSO", plotresult[[10]])
              saveWorkbook(wb,paste(dir,"/",i,"_result for plot.csv",sep=""), overwrite = TRUE)
            },silent=FALSE)
          }
          
          incProgress(1/(max(trait)-min(trait)+1), detail = paste("Doing part", i))
          Sys.sleep(0.1)
        }
      })
      
    }) 
    
    output$result<-renderUI(MR())
    
    ####################Draw plot###############################################
    plotfu<-reactive({
      Plot2<-function(fileplot=NULL,color1=NULL,color2=NULL,CriValue=NULL,p_stand=NULL,method=NULL,type=NULL){
        Manhattan<-function(fileplot,color1,color2,CriValue,method){
          parms<-as.data.frame(read.xlsx(fileplot,sheet = paste("Manhattan ",method,sep="")))
          mannewp<-as.numeric(parms[1,5])
          updateTextInput(session, "mancrimr", value = sprintf("%.6s",-log10(mannewp)))
          svgwline<-CriValue
          standline<-svgwline
          manhattan(parms,chr = "Chromosome",bp ="BPnumber",p ="P-value",snp="SNPname",col=c(color1,color2),suggestiveline=FALSE,genomewideline = standline)
        }
        QQplot1<-function(fileplot,p_stand,color1,color2,method){
          p_value<-as.matrix(read.xlsx(fileplot,sheet=paste("qq ",method,sep="")))
          pvalue<-matrix(p_value,,1)
          observed<-sort(pvalue[,1])
          observed<-observed/2
          observed<-observed[which(observed!=0)]
          newobserved<-observed[which(observed<(p_stand/2))]
          lobs<--(log10(newobserved))
          expected<-c(1:length(newobserved))
          lexp<--(log10(expected/(length(pvalue)+1)))
          plot(lexp,lobs,xlim=c(0,max(lexp)),ylim=c(0,max(lobs)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color2)
          abline(0,1,col=color1)
        }
        QQplot2<-function(fileplot,color1,color2,method){
          ress1<-as.data.frame(read.xlsx(fileplot,sheet = paste("qq ",method,sep="")))
          pvalue<-as.matrix(ress1)
          ps<-pvalue[,1]
          obs.x<-sort(ps)
          newobs.x<-obs.x[obs.x<1]
          n<-length(newobs.x)
          es<-(1:n)/(n+1)
          x<--log10(es)
          y<--log10(newobs.x)
          y<-y-0.3
          plot(x,y,xlim=c(0.3,max(x)),ylim=c(0.3,max(y)),xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'),col=color2)
          abline(0,1,col=color1)
        }
        LOD<-function(fileplot=NULL,color1,method=NULL){
          data<-as.data.frame(read.xlsx(fileplot,sheet=paste("Plot ",method,sep=""),colNames=TRUE))
          data<-sapply(data,as.numeric)
          gen<-data[,1:2]
          resulty<-data[,3:5]
          resultk<-which(is.na(data)==TRUE,arr.ind = TRUE)
          resultq<-resulty[1:(resultk[1]-1),] 
          if(nrow(resultq)>1){
            result<-resultq
          }else{
            result<-t(as.matrix(resultq))
          }
          galaxyy<-as.matrix(result)
          chr_pos <- gen[,1:2]
          chr_num <- length(unique(chr_pos[,1]))
          chr <- matrix(0,chr_num,1)
          pos <- matrix(0,chr_num,1)
          for(i in 1:chr_num)
          {
            temp <- numeric()
            temp <- length(which(chr_pos[,1]==i))
            if(i==1)
            {
              pos[i] <- temp
              chr[i] <- chr_pos[pos[i],2]
            }else{
              pos[i] <- pos[i-1] + temp
              chr[i] <- chr_pos[pos[i],2]
            }
          }
          pos_acc <- matrix(0,chr_num,1)
          for(i in 1:chr_num)
          {
            if(i==1){
              pos_acc[i] <- chr[i]
            }else{
              pos_acc[i] <- pos_acc[i-1] + chr[i]
            }
          }
          newres_pos <- galaxyy[,2]
          res_sumpos <- pos_acc[galaxyy[which(galaxyy[,1]>1),1]-1] + galaxyy[which(galaxyy[,1]>1),2] 
          newres_pos[which(galaxyy[,1]>1)] <- res_sumpos 
          pospic<-c(newres_pos)
          lodpic<-c(galaxyy[,3])  
          mm<-round(max(pospic)/4000)
          mm<-as.numeric(format(mm,digits = 1,scientific = TRUE))
          pospicx<-pospic/mm
          if(pospicx[1]<20){
            pospicx[1]<-pospicx[1]+20
          }
          pos_acc1<-pos_acc/mm
          resdf1 <- data.frame(pospicx,lodpic)
          pp <- ggplot(data=resdf1, aes(x=pospicx, y=lodpic)) +
            geom_bar(stat="identity", width=0.5, fill="white", linetype="solid",color=color1)
          
          pp <- pp + geom_vline(xintercept=c(0,pos_acc1),linetype="dashed",alpha=0.2)
          pp <- pp  + scale_x_continuous(expand=c(0,0),limits=c(0,(pos_acc1[dim(pos_acc1)[1]]+100))) +
            scale_y_continuous(expand=c(0,0))
          pp <- pp + xlab(paste("Genome position (",mm,"bp)",sep = "")) + ylab("LOD score") + ggtitle("") + theme_classic()
          pp <- pp + theme(axis.title.y = element_text( vjust = 2,hjust=0.5,size = 14),
                           axis.title.x = element_text(vjust = -0.5,hjust=0.5,size = 14))
          pp <- pp + theme(panel.background = element_rect(fill = "white"))
          pp <- pp + theme(text=element_text(family="mono"))
          pp <- pp + theme(axis.line.y = element_line(colour = "black", linetype = "solid"),
                           axis.line.x = element_line(colour = "black", linetype = "solid"))
          print(pp)
        }
        if(type=="Manhattan"){
          Manhattan(fileplot,color1,color2,CriValue,method)
        }else if(type=="qq"){
          if(method=="FASTmrEMMA"){
            QQplot2(fileplot,color1,color2,method)  
          }else{
            QQplot1(fileplot,p_stand,color1,color2,method) 
          }
        }else if(type=="LOD"){
          LOD(fileplot,color1,method)
        }
      }
    })
    
    output$mplmr<-renderPlot({
      plotfunction<-plotfu()
      req(input$filemanmr)
      plotfunction(fileplot=input$filemanmr$datapath,color1=input$manchcolour1mr,color2=input$manchcolour2mr,CriValue=as.numeric(input$mancrimr),p_stand=NULL,method=input$plotmethodman,type="Manhattan")
      
    })
    
    
    output$qplmr<-renderPlot({
      plotfunction<-plotfu()
      req(input$fileqqmr)
      plotfunction(fileplot=input$fileqqmr$datapath,color1=input$qqchcolour1mr,color2=input$qqchcolour2mr,CriValue=NULL,p_stand=0.95,method=input$plotmethodqq,type="qq")
      
    })
    
    
    output$qplPLA<-renderPlot({
      plotfunction<-plotfu()
      req(input$fileppPLA)
      plotfunction(fileplot=input$fileppPLA$datapath,color1=input$ppchcolour1PLA,color2=NULL,CriValue=NULL,p_stand=NULL,method=input$plotmethodlod,type="LOD")
      
    })
    
    output$downloadmanplot <- downloadHandler(
      filename = function() {
        paste("Manhattan", sep = ".", switch(
          input$plmanformat, "*.png"=".png", "*.tiff"=".tiff", "*.jpeg"=".jpeg","*.pdf"=".pdf"
        ))
      },
      content = function(file) {
        plotfunction<-plotfu()
        req(input$filemanmr)
        if(input$plmanformat=="*.png"){
          png(file,width=as.numeric(input$manwidthmr), height=as.numeric(input$manheimr), units= "px", pointsize =as.numeric(input$manwordremr),res=as.numeric(input$manfigureremr))
        }else if(input$plmanformat=="*.tiff"){
          tiff(file,width=as.numeric(input$manwidthmr), height=as.numeric(input$manheimr), units= "px", pointsize =as.numeric(input$manwordremr),res=as.numeric(input$manfigureremr))
        }else if(input$plmanformat=="*.jpeg"){
          jpeg(file,width=as.numeric(input$manwidthmr), height=as.numeric(input$manheimr), units= "px", pointsize =as.numeric(input$manwordremr),res=as.numeric(input$manfigureremr))
        }else if(input$plmanformat=="*.pdf"){
          pdf(file,width=10)
        }
        plotfunction(fileplot=input$filemanmr$datapath,color1=input$manchcolour1mr,color2=input$manchcolour2mr,CriValue=as.numeric(input$mancrimr),p_stand=NULL,method=input$plotmethodman,type="Manhattan")
        dev.off()
      })
    
    
    output$downloadqqplot <- downloadHandler(
      filename = function() {
        paste("QQ", sep = ".", switch(
          input$plqqformat, "*.png"=".png", "*.tiff"=".tiff", "*.jpeg"=".jpeg","*.pdf"=".pdf"
        ))
      },
      content = function(file) {
        plotfunction<-plotfu()
        req(input$fileqqmr)
        if(input$plqqformat=="*.png"){
          png(file,width=as.numeric(input$qqwidthmr), height=as.numeric(input$qqheimr), units= "px", pointsize =as.numeric(input$qqwordremr),res=as.numeric(input$qqfigureremr))
        }else if(input$plqqformat=="*.tiff"){
          tiff(file,width=as.numeric(input$qqwidthmr), height=as.numeric(input$qqheimr), units= "px", pointsize =as.numeric(input$qqwordremr),res=as.numeric(input$qqfigureremr))
        }else if(input$plqqformat=="*.jpeg"){
          jpeg(file,width=as.numeric(input$qqwidthmr), height=as.numeric(input$qqheimr), units= "px", pointsize =as.numeric(input$qqwordremr),res=as.numeric(input$qqfigureremr))
        }else if(input$plqqformat=="*.pdf"){
          pdf(file,width=10)
        }
        
        plotfunction(fileplot=input$fileqqmr$datapath,color1=input$qqchcolour1mr,color2=input$qqchcolour2mr,CriValue=NULL,p_stand=0.95,method=input$plotmethodqq,type="qq")
        
        dev.off()
      })
    
    
    
    output$downloadppplotPLA <- downloadHandler(
      filename = function() {
        paste("LOD", sep = ".", switch(
          input$plppformatPLA, "*.png"=".png", "*.tiff"=".tiff", "*.jpeg"=".jpeg","*.pdf"=".pdf"
        ))
      },
      content = function(file) {
        plotfunction<-plotfu()
        req(input$fileppPLA)
        if(input$plppformatPLA=="*.png"){
          png(file,width=as.numeric(input$ppwidthPLA), height=as.numeric(input$ppheiPLA), units= "px", pointsize =as.numeric(input$ppwordrePLA),res=as.numeric(input$ppfigurerePLA))
        }else if(input$plppformatPLA=="*.tiff"){
          tiff(file,width=as.numeric(input$ppwidthPLA), height=as.numeric(input$ppheiPLA), units= "px", pointsize =as.numeric(input$ppwordrePLA),res=as.numeric(input$ppfigurerePLA))
        }else if(input$plppformatPLA=="*.jpeg"){
          jpeg(file,width=as.numeric(input$ppwidthPLA), height=as.numeric(input$ppheiPLA), units= "px", pointsize =as.numeric(input$ppwordrePLA),res=as.numeric(input$ppfigurerePLA))
        }else if(input$plppformatPLA=="*.pdf"){
          pdf(file,width=10)
        }
        plotfunction(fileplot=input$fileppPLA$datapath,color1=input$ppchcolour1PLA,color2=NULL,CriValue=NULL,p_stand=NULL,method=input$plotmethodlod,type="LOD")
        dev.off()
      })
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