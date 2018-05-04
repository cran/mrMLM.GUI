server<-(function(input, output,session) {
  
  options(shiny.maxRequestSize=-1)
  
  output$general_ui <- renderUI({
    
    if (input$layout == "mrMLM") {
      
      tabPanel("",
               
               h2("Multi-locus GWAS methods",align="center"),
               br(),
               column(10,
                      h4("1. Zhang YM, Mao Y, Xie C, Smith H, Luo L, Xu S*. 2005. Mapping quantitative trait loci using 
                         naturally occurring genetic variance among commercial inbred lines of maize (Zea mays L.).
                         Genetics 169: 2267-2275"),
                      
                      h4("2. Wang SB, Feng JY, Ren WL, Huang B, Zhou L, Wen YJ, Zhang J, Jim M Dunwell, Xu S*, Zhang YM*.
                         2016. Improving power and accuracy of genome-wide association studies via a multi-locus mixed
                         linear model methodology. Scientific Reports 6:19444. (mrMLM)"),
                      
                      h4("3. Tamba CL. 2017. A fast mrMLM algorithm improves statistical power, accuracy and computational
                         efficiency of multi-locus genome-wide association studies. Nanjing Agricultural University. Ph D
                         Dissertation.(FASTmrMLM)"),
                      
                      h4("4. Zhang J, Feng JY, Ni YL, Wen YJ, Niu Y, Tamba CL, Yue C, Song QJ, Zhang YM*. pLARmEB: 
                         integration of least angle regression with empirical Bayes for multi-locus genome-wide association
                         studies. Heredity 118(6): 517-524.(pLARmEB)"),
                      
                      h4("5. Ren WL, Wen YJ, Jim M Dunwell, Zhang YM*. 2017. pKWmEB: integration of Kruskal-Wallis test with
                         empirical Bayes under polygenic background control for multi-locus genome-wide association study. 
                         Heredity, online (2017-12-13), doi:10.1038/s41437-017-0007-4.(pKWmEB)"),             
                      
                      h4("6. Wen YJ, Zhang H, Ni YL, Huang B, Zhang J, Feng JY, Wang SB, Jim M Dunwell, Zhang YM*, Wu R*. 2017.
                         Methodological implementation of mixed linear models in multi-locus genome-wide association studies. 
                         Briefings in Bioinformatics, bbw145, https://doi.org/10.1093/bib/bbw145.(FASTmrEMMA)"),
                      
                      h4("7. Tamba CL, Ni YL, Zhang YM*. Iterative sure independence screening EM-Bayesian LASSO algorithm for
                         multi-locus genome-wide association studies. PLoS Computational Biology 2017,13(1):e1005357.
                         (ISIS EM-BLASSO)"),
                      
                      br(),
                      br(),
                      
                      h4("Authors: Zhang Ya-Wen,Li Pei, Ren Wen-Long, Ni Yuan-Li, Zhang Yuan-Ming"),
                      
                      h4("Maintainer: Zhang Yuan-Ming (soyzhang at mail.hzau.edu.cn)"), 
                      
                      h4("mrMLM version 3.1, Realeased March 2018"),
                      
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
                            
                            radioButtons("radiomrformat", "Select dataset format",
                                         choices = list("mrMLM numeric format", "mrMLM character format","Hapmap (TASSEL) format"),selected ="mrMLM numeric format"),
                            offset = 1
                     ),
                     column(3,
                            br(),
                            fileInput("filegenmr", "Input Genotypic file",multiple = TRUE,accept = ".csv")
                     ),
                     
                     column(3,
                            br(),
                            radioButtons("dispgenmr", "Display genotype", choices = c(Head = "head",All = "all"),selected = "head"),
                            offset = 1
                     ),
                     column(12,
                            dataTableOutput("contentsrawgenmr"),
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
                            fileInput("filephenmr", "Input Phenotypic file",multiple = TRUE,accept = ".csv"),
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
                            radioButtons("ckmr", "Input Kinship?",
                                         choices = list("Input Kinship (K) matrix file"="inmr", "Calculate Kinship (K) matrix by this software"="cmr"),selected = "inmr"), 
                            offset = 1
                     ),
                     column(3,
                            br(),
                            fileInput("filekmr", "Input Kinship (K)",multiple = TRUE,accept = ".csv"),
                            helpText("Note:Please select the 1st option if no. of markers is more than 50,000."),
                            offset = 1
                            
                     ),
                     
                     column(3,
                            br(),
                            radioButtons("dispkinmr", "Display kinship", choices = c(Head = "head",All = "all"),selected = "head"),
                            offset = 1
                     ),
                     
                     column(12,
                            br(),
                            dataTableOutput("kinshipmr"),
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
                     
                     column(3,
                            br(),
                            radioButtons("instrmr", "Include population struction(Q) matrix?", choices = c("Not included in the model"="notinc" ,"Included"="inc" ),selected = "notinc"), 
                            offset = 1
                     ),
                     column(3,
                            br(),
                            fileInput("filestrmr", "Input population structure",multiple = TRUE,accept = ".csv"),
                            offset = 1
                            
                     ),
                     
                     
                     column(3,
                            br(),
                            radioButtons("dispstrmr", "Display population structure", choices = c(Head = "head",All = "all"),selected = "head"),
                            offset = 1
                     ),
                     
                     
                     column(12,
                            br(),
                            dataTableOutput("psmr"),
                            offset = 1
                     )
                   )
                   
          ),
          
          tabPanel("Method select & Parameter settings",
                   
                   fluidRow(  
                     column(12,
                            checkboxGroupInput("Method", "Method select",
                                               c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),selected="mrMLM",inline=TRUE),
                            
                            offset = 1
                     ),
                     column(3,
                            br(),
                            textInput("lodmr", label = "Critical LOD score (All method)",value = "3"),
                            offset = 1
                            
                     ),
                     column(3,
                            textInput("scgmr", label = "Search radius of candidate gene (kb)(mrMLM & FASTmrmLM):",value = "20")
                            
                     ),
                     
                     column(3,
                            textInput("scgPLA", label = "No. of potentially associated variables selected by LARS (pLARmEB):",value = "50")
                     ),
                     
                     column(3,
                            radioButtons("chdofunFME", "Likelihood Function(FASTmrEMMA)", choices = c("REML","ML"),inline=TRUE),
                            offset = 1
                     ),
                     
                     column(3,
                            radioButtons("BootPLA", "Bootstrap (pLARmEB)", choices = c("TRUE","FALSE"),inline=TRUE,selected="FALSE")
                     ),
                     
                     column(3,
                            radioButtons("Drawall", "Draw plot or not (All method)", choices = c("TRUE","FALSE"),inline=TRUE,selected="FALSE")
                            
                     ),
                     
                     column(3,
                            radioButtons("Resolution", "Plot resolution (All method)", choices = c("High", "Low"),inline = TRUE),
                            offset = 1
                     ),
                     
                     
                     column(3,
                            radioButtons("Plotformat", "Plot format (All method)", choices = c("*.png", "*.tiff", "*.jpeg","*.pdf"),inline = TRUE)
                     ),
                     
                     column(3,
                            textInput("SavePath", "Save path",value = "C:/Users/Administrator/Desktop/")
                           ),
                     
                     
                     column(3,
                            h4("Please select trait ID") ,
                            textInput("Trait1", "From", value="1"),
                            offset = 1
                     ),
                     
                     column(3,
                            
                            br(),
                            br(),
                            
                            textInput("Trait2", "To", value="1")
                     ),
                     
                     
                     column(12, 
                            
                            br(),
                            br(),
                            br(),
                            actionButton("runmr", label = "Run",width=280, icon("paper-plane"), 
                                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                            
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
                            radioButtons("Precisionman", "Select resolution of plot", c("General precision", "High precision", "Set by yourself"),inline = TRUE)
                     ),
                     column(3,
                            textInput("manwidthmr", label = "Figure Width (px):",value = "960")
                            
                     ),
                     column(3,
                            textInput("manheimr", label = "Figure height (px)",value = "600")
                     ),
                     column(3,
                            textInput("manwordremr", label = "Word resolution (1/72 inch, ppi):",value = "20")
                     ),
                     
                     column(3,
                            textInput("manfigureremr", label = "Figure resolution (ppi):",value = "72")
                     ),
                     column(3,
                            selectInput(inputId = "manchcolour1mr",label = "Chromosome color (odd):",
                                        choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "blue")
                     ),
                     column(3,
                            selectInput(inputId = "manchcolour2mr",label = "Chromosome color (even):",
                                        choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "red")
                     ),
                     
                     column(3,
                            textInput("mancrimr", label = "Critical value for Manhattan Plot:",value = "3")
                     ),
                     
                     column(3,
                            selectInput(inputId = "plotmethodman",label = "Method select",
                                        choices = c("mrMLM","FASTmrMLM","FASTmrEMMA","pKWmEB"),selected = "mrMLM")
                     ),
                     
                     column(12,
                            fileInput("filemanmr", "Input file to draw Manhattan plot",multiple = TRUE,accept = ".csv")
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
                            radioButtons("Precisionqq", "Select resolution of plot", c("General precision", "High precision", "Set by yourself"),inline = TRUE) 
                     ),
                     
                     column(3,
                            textInput("qqwidthmr", label = "Figure Width (px):",value = "960")
                     ),
                     column(3,
                            textInput("qqheimr", label = "Figure height (px)",value = "600")
                     ),
                     column(3,
                            textInput("qqwordremr", label = "Word resolution (1/72 inch, ppi):",value = "20")
                     ),
                     
                     column(3,
                            textInput("qqfigureremr", label = "Figure resolution (ppi):",value = "72")
                     ),
                     column(3,
                            selectInput(inputId = "qqchcolour1mr",label = "Point color:",
                                        choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "red")
                     ),
                     column(3,
                            selectInput(inputId = "qqchcolour2mr",label = "Line color:",
                                        choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "black")
                     ),
                     
                     column(3,
                            selectInput(inputId = "plotmethodqq",label = "Method select",
                                        choices = c("mrMLM","FASTmrMLM","FASTmrEMMA","pKWmEB"),selected = "mrMLM")
                     ),
                     
                     column(12,
                            fileInput("fileqqmr", "Input file to draw QQ plot",multiple = TRUE,accept = ".csv")
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
                            radioButtons("Precisionlod", "Select resolution of plot", c("General precision", "High precision", "Set by yourself"),inline = TRUE) 
                     ),
                     
                     column(4,
                            textInput("ppwidthPLA", label = "Figure Width (px):",value = "960")
                     ),
                     column(4,
                            textInput("ppheiPLA", label = "Figure height (px)",value = "240")
                     ),
                     column(4,
                            textInput("ppwordrePLA", label = "Word resolution (1/72 inch, ppi):",value = "12")
                     ),
                     
                     column(4,
                            textInput("ppfigurerePLA", label = "Figure resolution (ppi):",value = "72")
                     ),
                     column(4,
                            selectInput(inputId = "ppchcolour1PLA",label = "Lod line color:",
                                        choices = c("blue","black","red","yellow","green","pink","purple","gray","brown"),selected = "red")
                     ),
                     
                     column(4,
                            selectInput(inputId = "plotmethodlod",label = "Method select",
                                        choices = c("pLARmEB","ISIS EM-BLASSO"),selected = "pLARmEB")
                     ),
                     
                     
                     column(12,
                            fileInput("fileppPLA", "Input file to draw plot",multiple = TRUE,accept = ".csv")
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
  
  
  manual<-eventReactive(input$manl,{
    
    RShowDoc("Instruction",package="mrMLM.GUI") 
    
  })
  
  output$manl<-renderUI(manual())
  
  
  
  
  
  
  upda1<-observeEvent(input$Precisionman,{
    
    if(input$Precisionman=="High precision"){
      widthman<-10000
      heightman<-6000
      wordman<-30
      figureman<-300
    }else if(input$Precisionman=="General precision"){
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
    
    if(input$Precisionqq=="High precision"){
      widthqq<-10000
      heightqq<-6000
      wordqq<-30
      figureqq<-300
    }else if(input$Precisionqq=="General precision"){
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
    if(input$Precisionlod=="High precision"){
      widthlod<-10000
      heightlod<-6000
      wordlod<-30
      figurelod<-600
    }else if(input$Precisionlod=="General precision"){
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
  
  
  output$contentsrawgenmr <- renderDataTable({
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
    kkRaw[1,2:nnkk]<-"  "
    kkRaw
    
  })
  
  output$kinshipmr<-renderDataTable({
    
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
  
  strcmrshow<-reactive({
    psmatrixRaw<-strcmr()
    nnpprow<-dim(psmatrixRaw)[1]
    nnppcol<-dim(psmatrixRaw)[2]
    psmatrixRaw[1,2:nnppcol]<-"  "
    psmatrixRaw
    
  })
  
  
  output$psmr<-renderDataTable({
    if(input$instrmr=="inc"){
      if(input$dispstrmr == "head"){
        return(head(strcmrshow()))
      }else{
        return(strcmrshow())
      }
      
    } 
    
  })
  
  observe({
    shinyjs::toggleState("filestrmr",input$instrmr!="notinc")
  })
  
  
  
  ReadData<-reactive({
    kkRaw<-NULL
    psmatrixRaw<-NULL
    genRaw<-genRawmr()
    pheRaw1q<- pheRawmr1()
    
    if(input$ckmr=="inmr"&&input$filekmr!=""){
      kkRaw<-kinmrshow()
    }
    if(input$instrmr=="inc"&&input$filestrmr !=""){
      psmatrixRaw<-strcmr()
    }
    phename<-as.matrix(pheRaw1q[1,2:ncol(pheRaw1q)])
    output<-list(genRaw=genRaw,pheRaw1q=pheRaw1q,kkRaw=kkRaw,psmatrixRaw=psmatrixRaw,phename=phename)
  })
  
  
  MR<-eventReactive(input$runmr,{
    
    dir<-input$SavePath
    setwd(dir)
    
    DoData<-function(genRaw=NULL,Genformat=NULL,pheRaw1q=NULL,kkRaw=NULL,psmatrixRaw=NULL,trait=NULL,type=NULL){
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
        genRaw[which(genRaw==0)]<-0.5
        genRaw[which(genRaw==-1)]<-0 
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
        gen<-as.matrix(needGen[-1,])
        gen<-matrix(as.numeric(gen),nrow=nrow(gen))
        ###########To show on the table ----newPhe
        pheRaw[1,2]<-"  "
        newPhe<-rbind(pheRaw[1,],newPhe)
        ###########To be computed ----phe
        phe<-as.matrix(newPhe[-1,-1])
        phe<-matrix(as.numeric(phe),nrow=nrow(phe))
        shownewGen<-newGen[-1,]
        colnames(shownewGen)<-newGen[1,]
        shownewGen<-as.data.frame(shownewGen)
        shownewPhe<-newPhe[-1,]
        colnames(shownewPhe)<-c(newPhe[1,1],"   ")
        shownewPhe<-as.data.frame(shownewPhe)
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
        ###########To be computed ----gen
        locChr <- as.numeric(which(newGen[1,]=="chrom"))
        locPos <- as.numeric(which(newGen[1,]=="pos"))
        needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
        needGen<-newGen[,needloc]
        gen<-as.matrix(needGen[-1,])
        gen<-matrix(as.numeric(gen),nrow=nrow(gen))
        ###########To show on the table ----newPhe
        pheRaw[1,2]<-"  "
        newPhe<-rbind(pheRaw[1,],newPhe)
        ###########To be computed ----phe
        phe<-as.matrix(newPhe[-1,-1])
        phe<-matrix(as.numeric(phe),nrow=nrow(phe))
        shownewGen<-newGen[-1,]
        colnames(shownewGen)<-newGen[1,]
        shownewGen<-as.data.frame(shownewGen)
        shownewPhe<-newPhe[-1,]
        colnames(shownewPhe)<-c(newPhe[1,1],"   ")
        shownewPhe<-as.data.frame(shownewPhe)
        
      }else if(inputform==3){
        ##########To find the same individual ID between genotype and phenotype
        nameGen<-as.matrix(genRaw[1,],1,)
        namePhe<-as.matrix(pheRaw[,1],,1)
        sameName<-intersect(nameGen,namePhe)
        ##########To find the location of the same name 
        locGen<-match(sameName,nameGen)
        locPhe<-match(sameName,namePhe)
        ##########Produce new genotype matrix and phenotype matrix
        hapName<-matrix(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode"),1,)
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
        ###########To be computed ----gen
        locChr<-as.numeric(which(newGen[1,]=="chrom"))
        locPos<-as.numeric(which(newGen[1,]=="pos"))
        needloc<-c(locChr,locPos,(nnhap+1):colnewGen)
        needGen<-newGen[,needloc]
        gen<-as.matrix(needGen[-1,])
        gen<-matrix(as.numeric(gen),nrow=nrow(gen))
        ###########To show on the table ----newPhe
        pheRaw[1,2]<-"  "
        newPhe<-rbind(pheRaw[1,],newPhe)
        ###########To be computed ----phe
        phe<-as.matrix(newPhe[-1,-1])
        phe<-matrix(as.numeric(phe),nrow=nrow(phe))
        shownewGen<-newGen[-1,]
        colnames(shownewGen)<-newGen[1,]
        shownewGen<-as.data.frame(shownewGen)
        shownewPhe<-newPhe[-1,]
        colnames(shownewPhe)<-c(newPhe[1,1],"   ")
        shownewPhe<-as.data.frame(shownewPhe)
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
        ##revised
        filtername<-as.vector(psmatrixRaw[2,2:nnppcol])
        selectpsmatrix<-matrix(as.numeric(psmatrixPre[locPop,-1]),nrow = length(locPop))
        psum<-apply(selectpsmatrix,1,sum)
        psum<-round(psum)
        sumps<-sum(psum)
        m<-dim(selectpsmatrix)[1]
        if(sumps>=m){
          combovalue<-"Q1"
          coldelet<-unlist(str_extract_all(combovalue,"[0-9]+"))
          coldelet<-as.numeric(coldelet)
          psmatrix<-as.matrix(selectpsmatrix[,-coldelet])
          psmatrixRaw<-as.matrix(psmatrixRaw[,-(coldelet+1)])
        }else{
          psmatrix<-selectpsmatrix
        }
      }
      
      doresult<-list(gen=gen,phe=phe,outATCG=outATCG,genRaw=genRaw,kk=kk,psmatrix=psmatrix)
      return(doresult)
    }
    
    
    InputDataF<-function(readraw,Genformat=NULL,method=NULL,trait=NULL){
      
      doMR<-NULL;doFME<-NULL
      
      if("mrMLM"%in%method){
        doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2)
      }
      
      if("FASTmrMLM"%in%method){
        doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2)
      }  
      
      if("FASTmrEMMA"%in%method){
        doFME<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=1)  
      }
      
      if("pLARmEB"%in%method){
        doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2) 
      }
      if("pKWmEB"%in%method){
        doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2) 
      }  
      
      if("ISIS EM-BLASSO"%in%method){
        doMR<-DoData(readraw$genRaw,Genformat,readraw$pheRaw1q,readraw$kkRaw,readraw$psmatrixRaw,trait,type=2) 
      }
      
      output<-list(doMR=doMR,doFME=doFME) 
      return(output)
      
    }
    
    
    svrad<-as.numeric(input$scgmr)
    svmlod<-as.numeric(input$lodmr)
    lars1<-as.numeric(input$scgPLA)
    Likelihood<-input$chdofunFME
    Bootstrap<-input$BootPLA 
    DrawPlot<-input$Drawall
    Plotformat<-input$Plotformat
    Resolution<-input$Resolution
    
    
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
    
    trait1<-as.numeric(input$Trait1);trait2<-as.numeric(input$Trait2)
    
    print("Running in progress, please be patient...")
    
    
    withProgress(message = 'Running in progress', value = 0, {
      
      
      for (i in trait1:trait2){
        
        InputData<-InputDataF(readraw,Genformat,input$Method,i)
        
        reMR<-NULL;reFMR<-NULL;reFME<-NULL;rePLA<-NULL;rePKW<-NULL;reISIS<-NULL
        re1MR<-NULL;re1FMR<-NULL;re1FME<-NULL;re1PLA<-NULL;re1PKW<-NULL;re1ISIS<-NULL
        remanMR<-NULL;reqqMR<-NULL;remanFMR<-NULL;reqqFMR<-NULL;remanFME<-NULL;reqqFME<-NULL;replPLA<-NULL;remanPKW<-NULL;reqqPKW<-NULL; replISIS<-NULL
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
                png(paste(i,"_mrMLM_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              }else if(Plotformat1=="*.tiff"){
                tiff(paste(i,"_mrMLM_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              }else if(Plotformat1=="*.jpeg"){
                jpeg(paste(i,"_mrMLM_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              }else if(Plotformat1=="*.pdf"){
                pdf(paste(i,"_mrMLM_Manhattan.pdf"),width=10)
              }
              
              Plot(plotresult=remanMR,color1="red",color2="blue",0.95,method="mrMLM",type="Manhattan")
              dev.off()
              
              if(Plotformat2=="*.png"){
                png(paste(i,"_mrMLM_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              }else if(Plotformat2=="*.tiff"){
                tiff(paste(i,"_mrMLM_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              }else if(Plotformat2=="*.jpeg"){
                jpeg(paste(i,"_mrMLM_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
              }else if(Plotformat2=="*.pdf"){
                pdf(paste(i,"_mrMLM_qq.pdf"),width=10)
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
                  png(paste(i,"_FASTmrMLM_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(i,"_FASTmrMLM_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(i,"_FASTmrMLM_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(i,"_FASTmrMLM_Manhattan.pdf"),width=10)
                }
                
                Plot(plotresult=remanFMR,color1="red",color2="blue",0.95,method="FASTmrMLM",type="Manhattan")
                dev.off()
                
                if(Plotformat2=="*.png"){
                  png(paste(i,"_FASTmrMLM_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.tiff"){
                  tiff(paste(i,"_FASTmrMLM_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.jpeg"){
                  jpeg(paste(i,"_FASTmrMLM_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.pdf"){
                  pdf(paste(i,"_FASTmrMLM_qq.pdf"),width=10)
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
                  png(paste(i,"_FASTmrEMMA_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(i,"_FASTmrEMMA_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(i,"_FASTmrEMMA_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(i,"_FASTmrEMMA_Manhattan.pdf"),width=10)
                }
                
                Plot(plotresult=remanFME,color1="red",color2="blue",0.95,method="FASTmrEMMA",type="Manhattan")
                dev.off()
                
                if(Plotformat2=="*.png"){
                  png(paste(i,"_FASTmrEMMA_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.tiff"){
                  tiff(paste(i,"_FASTmrEMMA_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.jpeg"){
                  jpeg(paste(i,"_FASTmrEMMA_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.pdf"){
                  pdf(paste(i,"_FASTmrEMMA_qq.pdf"),width=10)
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
              }
              replPLA<-outPLA$plot
              
              if(DrawPlot==TRUE){
                if(Resolution=="Low"){
                  manwidth<-960;manhei<-240;manwordre<-12;manfigurere<-72
                }else if(Resolution=="High"){
                  manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-600
                }
                if(Plotformat1=="*.png"){
                  png(paste(i,"_pLARmEB_LOD.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(i,"_pLARmEB_LOD.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(i,"_pLARmEB_LOD.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(i,"_pLARmEB_LOD.pdf"),width=12)
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
                  png(paste(i,"_pKWmEB_Manhattan.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(i,"_pKWmEB_Manhattan.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(i,"_pKWmEB_Manhattan.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(i,"_pKWmEB_Manhattan.pdf"),width=10)
                }
                
                Plot(plotresult=remanPKW,color1="red",color2="blue",0.95,method="pKWmEB",type="Manhattan")
                dev.off()
                
                if(Plotformat2=="*.png"){
                  png(paste(i,"_pKWmEB_qq.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.tiff"){
                  tiff(paste(i,"_pKWmEB_qq.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.jpeg"){
                  jpeg(paste(i,"_pKWmEB_qq.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat2=="*.pdf"){
                  pdf(paste(i,"_pKWmEB_qq.pdf"),width=10)
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
              }
              replISIS<-outISIS$plot
              
              if(DrawPlot==TRUE){
                if(Resolution=="Low"){
                  manwidth<-960;manhei<-240;manwordre<-12;manfigurere<-72
                }else if(Resolution=="High"){
                  manwidth<-10000;manhei<-6000;manwordre<-30;manfigurere<-600
                }
                if(Plotformat1=="*.png"){
                  png(paste(i,"_ISIS EM-BLASSO_LOD.png"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.tiff"){
                  tiff(paste(i,"_ISIS EM-BLASSO_LOD.tiff"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.jpeg"){
                  jpeg(paste(i,"_ISIS EM-BLASSO_LOD.jpeg"),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
                }else if(Plotformat1=="*.pdf"){
                  pdf(paste(i,"_ISIS EM-BLASSO_LOD.pdf"),width=12)
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
            
            write.table(output,paste(i,"_Final result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
            
            write.table(output1,paste(i,"_intermediate result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
            
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
            
            saveWorkbook(wb,paste(i,"_resultforplot.xlsx",sep=""), overwrite = TRUE)
          },silent=FALSE)
        }
        incProgress(1/(trait2-trait1+1), detail = paste("Doing part", i))
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