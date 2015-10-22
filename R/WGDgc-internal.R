.calcProb <- function (logLamlogMu, ti, mStart, nEnd)
{ # calculate probability from mStart genes in a parent to nEnd genes in a child in time=ti	
  if (mStart==0 && nEnd==0) {return(1)}
  if ( (length(logLamlogMu)==1) || (logLamlogMu[1] == logLamlogMu[2]) )
    return( .calcProbSamePara(logLamlogMu[1], ti, mStart, nEnd) )
  logLam = logLamlogMu[1]; logMu = logLamlogMu[2];
  lmdiff = exp(logLam) - exp(logMu)	# lam - mu
  mlratio = exp(logMu) / exp(logLam)	# mu/lam
  beta = (exp(lmdiff *ti) -1) / (exp(lmdiff *ti) -mlratio)
  alpha = beta * mlratio
  minMN = min(mStart,nEnd)
  value = choose(mStart+nEnd-1, mStart-1) *alpha^mStart *beta^nEnd 	# j=0 in the formula
  Pmn = value
  ab = (1-alpha-beta)/(alpha*beta)	
  if (minMN >=1)
    for (j in 1:minMN) {
      value = value *(mStart-j+1)/j *(nEnd-j+1)/(mStart+nEnd-j) *ab
      Pmn = Pmn + value
    }
  return (Pmn)
}
.calcProbOneBranch <-
  function(branchNode, logLamlogMu, nPos, nFamily, nLeaf, 
           geneCountData, wgdTab, phyloMat, MatDoomed, MDcol, Mat)
{# branchNode is a vector of nodes on a branch in correct order as on the tree
 #  ex:(6,15,14,13,12,7) denotes 6 is parent of 15, 15 is parent of 14 ... 12 is parent of 7	 
 # branchNode[1] = ancestor of this branch, it is internal-nonsingleton node
 # branchNode[length(branchNode)] = last decendants of this branch, it is either 
 #  a leaf or internal-nonsingleton node
	 
 # if length(branchNode) >2, then this branch has some singletons in the middle
 # MDcol is the column of MatDoomed, either 2 or 3
 # this function returns an nPos x nFamily matrix of prob of data below ancestor node of this branch
  for(i in (length(branchNode) -1) :1)
    { # assume the prob of data below the last decendant is available at this point
      # so calculate from the node next to the last up to the oldest ancestor		
      parent = branchNode[i]
      child = branchNode[i+1]
      time = phyloMat$Time[phyloMat$Child==child]
      prob = matrix(0, nrow=nPos, ncol=nFamily)		
      if (child <= nLeaf) { # this child is a leaf
        species = phyloMat$Species [phyloMat$Child == child]
        prob = .getPt(logLamlogMu=logLamlogMu, ti=time, nPos=nPos, nFamily=nFamily, isChildLeaf=T,
          nLeafCount=geneCountData[ ,species])
        MatDoomed[parent,MDcol]=prob[2,nFamily]
        if (length(branchNode)>2) { #the parent is not an internal-nonSingleton node			
          MatDoomed[parent,ifelse(MDcol==2,3,2)]=1
          MatDoomed[parent,1]=MatDoomed[parent,MDcol] 			
        }
        # getPt = function (logLamlogMu, ti, nPos, nFamily=NULL, isChildLeaf=F, nLeafCount=NULL, 
	#   isWgdNode=F, wgdLR=NULL) {					
      } else { # child is an internal node
        if (time == 0) { # this parent is the node before wgd
          wgdLR = wgdTab$LossRate[wgdTab$wgdNode==parent]			
          Pt = .getPt (logLamlogMu=logLamlogMu, ti=time, nPos=nPos, isWgdNode=T, wgdLR=wgdLR)
          MatDoomed[parent,ifelse(MDcol==2,3,2)]=1
          MatDoomed[parent,MDcol]= wgdLR*MatDoomed[child,1] +(1-wgdLR)*MatDoomed[child,1]^2
          MatDoomed[parent,1]=MatDoomed[parent,MDcol]
        } else { # parent is not node before WGD and child is not a leaf
          Pt = .getPt (logLamlogMu=logLamlogMu, ti=time, nPos=nPos)
          MatDoomed[parent,MDcol]=Pt[2,1]			
          for (j in 2:nPos)
            MatDoomed[parent,MDcol]=MatDoomed[parent,MDcol]+ Pt[2,j] * MatDoomed[child,1]^(j-1)
          if(i!=1){
            MatDoomed[parent,ifelse(MDcol==2,3,2)]=1
            MatDoomed[parent,1]=MatDoomed[parent,MDcol]
          }
        }			
        for (k in 1:nPos)
          prob[k, ] = Pt[k, ] %*% Mat[ ,child-nLeaf, ]
      }		
      Mat[ ,parent-nLeaf, ] = prob
    }
  return(list(prob=prob, MatDoomed=MatDoomed, Mat=Mat))
}
.calcProbSamePara <- function (logLam, ti, mStart, nEnd)
{ # calculate probability from mStart genes in a parent to nEnd genes in a child in time=ti
  # here lamda = birth rate = death rate	
  if (mStart==0 && nEnd==0)  return(1)
  alpha = (exp(logLam) *ti) / (1+ exp(logLam)*ti)
  minMN = min(mStart, nEnd)
  value = choose(mStart+nEnd-1, mStart-1) *alpha^(mStart+nEnd)	# j=0 in the formula
  Pmn = value
  ab = (1-2*alpha) / alpha^2	
  if (minMN >= 1) 
    for (j in 1:minMN){
      value = value *(mStart-j+1)/j *(nEnd-j+1)/(mStart+nEnd-j) *ab
      Pmn = Pmn + value
    }
  return (Pmn)
}

.getPt <-
function (logLamlogMu, ti, nPos, nFamily=NULL, isChildLeaf=F, nLeafCount=NULL, isWgdNode=F, wgdLR=NULL) {
	# nLeafCount is 1 column in the geneCountData, it is a vector of genecounts of 1 species in all families
	# nPos = max posible number of genes at each internal node, for ex: nPos=101 (0 to 100 genes) 
	#   recommend choosing nPos=2*max(geneCountData)
		
	# this function return Pt which is a probability matrix
	
	# if isChildLeaf=T, Pt is a nPos x nFamily matrix in which each entry Pt[i,j]=prob from i-1 genes 
	#   to nLeafCount[j] genes in family j 
	# if isChildLeaf=F, Pt is nPos x nPos matrix, each entry Pt[i,j] = prob from i-1 genes to j-1 genes
	#   (based on birth-death process)
	
	# if isWgdNode=T, Pt is nPos x nPos matrix, each entry Pt[i,j] = prob from i-1 genes to j-1 genes
	#   (based on binomial distribution)
	# isWgdNode = whether or not this is a node before WGD event
	# wgdLR = wgd loss rate
  if (isWgdNode) {
    Pt = matrix (0, nrow=nPos, ncol=nPos) 
    for (i in 0:(nPos-1)) {
      N = ifelse( 2*i > (nPos-1), nPos-1, 2*i)			
      for (j in i:N) {
        Pt[i+1,j+1] = choose(i, j-i) *(1-wgdLR)^(j-i) *wgdLR^(2*i -j)
        # binomial dist of getting j genes (after wgd) from i genes (before wgd),
	# j is between i:2i or i:(nPos-1)
      }
    }
    return (Pt)
  }		
  if (isChildLeaf) {
    Pt = matrix(0, nrow=nPos, ncol=nFamily)
    for (i in 1:nPos) {
      for (j in 1:nFamily) {
        Pt[i,j] = .calcProb (logLamlogMu, ti, i-1, nLeafCount[j])
      }
    }
  } else {
    Pt = matrix (0, nPos, nPos)
    for (i in 1:nPos) {
      for (j in 1:nPos) {
        Pt[i,j] = .calcProb (logLamlogMu, ti, i-1, j-1)
      }
    }
  }
  # calcProb = function (logLamlogMu, ti, mStart, nEnd) {}
  return (Pt)
}


.checkClades <- function(phyloMat, geneCountData, nLeaf)
{ # check if all families have at least 1 gene copy in each clade
  r = which(phyloMat$Species == "Root")
  descRoot=which(phyloMat$Parent==r)
  tableClade<-array (0, nLeaf)
  for (Leafnode in 1 : nLeaf) {
    node<-Leafnode
    if (node==descRoot[1])
      tableClade[Leafnode] =1 
    else if (node==descRoot[2])
      tableClade[Leafnode] = -1
    while (node!=descRoot[1] & node!=descRoot[2]) { 
      indnode<-which(phyloMat$Child==node)
      node<-phyloMat$Parent[indnode]
      if (node==descRoot[1])
        tableClade[Leafnode] =1 
      else if (node==descRoot[2])
        tableClade[Leafnode] = -1
    }
  }
  codeClade=-1
  ## we begin with the first clade
  for (clade in 1:2){ # loop on the different clades
    indClade=which(tableClade==codeClade)
    SpecName<-phyloMat$Species[indClade]
    ind<-which(colnames(geneCountData)%in%SpecName)
    geneCountDataClade <- geneCountData[, c(ind)]
    if (length(ind)>1) 
      geneNumberInEachFamilyClade<-rowSums(geneCountDataClade, dims = 1)
    else
      geneNumberInEachFamilyClade<-geneCountDataClade
    indEmptyFamily=which(geneNumberInEachFamilyClade==0)
    if (length(indEmptyFamily)>0){
      print("The following families are responsible for the ERROR")
      print(indEmptyFamily)
      stop("ERROR about the conditioning, we can not use the type oneInBothClades 
		     	since at least one family has 0 gene copies in one of the 2 clades")
    } else {
      codeClade=1 # we have to check the second clade
    }
  }
}


.Random.seed <-
c(403L, 20L, -1637841155L, -723990537L, -742219490L, -1608373288L, 
-15061701L, 1375759065L, -2141665688L, -1658284138L, -1663114479L, 
-219566909L, 2078508242L, -2076888732L, 1944560807L, -1190172163L, 
534250484L, 993909978L, -746467803L, 275661631L, -745991994L, 
600403280L, -567683725L, 1638842257L, 776655296L, -270216290L, 
585289065L, 1712493979L, -2082539350L, -1315576116L, -335607921L, 
-2037824571L, -171153540L, -174827566L, 482113357L, -304077305L, 
1089320430L, 379076808L, 1297846219L, -1248744919L, 929781880L, 
-2112727642L, -2028003391L, 442631251L, -157947134L, -1629049292L, 
1179823031L, 306871213L, 171897028L, -409533206L, 683199317L, 
-694944529L, -2140646666L, -1295573472L, 389538787L, -1510463359L, 
-1909164304L, 1857672590L, -706946311L, -1491887989L, -1430839110L, 
-438520900L, -252843777L, 889404501L, -701249684L, 1479145410L, 
805404701L, -267140969L, -2089828866L, -765620808L, -547798629L, 
115083257L, -574566584L, 1830017654L, 1916867505L, 2078215651L, 
1670861490L, -1881821436L, -331115449L, -1532583139L, 1245690260L, 
1737194810L, 1578918277L, 764555871L, -2142022426L, 908896112L, 
-1939904365L, 1585095985L, 1613067232L, 720281726L, 1038646025L, 
-2134139589L, -874204918L, -656908052L, 1337347439L, 2049774245L, 
-679695268L, -879612878L, -358061267L, 177975783L, 820038926L, 
2134401768L, -1727976213L, -686004855L, 268809496L, 663320646L, 
536404129L, -2028525773L, -1255229790L, -2013395564L, -1508946921L, 
-2112758771L, 1841821988L, 36427018L, 1009235061L, 760176591L, 
-1916534570L, -599474048L, 1090584259L, 589187169L, -982124336L, 
-412372434L, -1814375463L, -1601905685L, -2142160294L, 1617382428L, 
1567463711L, -1538651787L, 712738828L, 1617383906L, -2032194883L, 
-1892326089L, 1912654302L, 983715608L, -1999174021L, 1254966297L, 
-156175192L, -1857368490L, -776639791L, -1706194301L, -1947385582L, 
-1031412572L, -745917337L, 1402191293L, 1811627828L, -159556454L, 
958374629L, -1628480385L, 1293604358L, -1809817200L, -135390157L, 
-1047415855L, -2106328704L, -1354181538L, -2142567511L, -1748149669L, 
1765621354L, 1378190092L, 1051223375L, -1781850619L, 601334972L, 
-697444974L, -697725299L, 96404679L, -1484004946L, -991900664L, 
1265565707L, 1826239721L, 1500402616L, -1347576730L, -963740415L, 
261954451L, 56081474L, 1345127412L, 822447607L, 1747309037L, 
-1321672060L, 23319722L, -931511787L, 196525743L, 1472134326L, 
315280096L, -425920989L, -1305282239L, 1376740400L, 1251638862L, 
-766645575L, 1927728843L, -386091142L, 1687043708L, -1436566849L, 
386456469L, 654056364L, 1871516034L, -336865443L, 377633367L, 
-1565675202L, 2134757752L, -1906365861L, -60891719L, -1337327608L, 
204620470L, -612738575L, 1017122851L, -1634652046L, -1595862332L, 
115606663L, -1749567395L, 2111921748L, 1569889402L, -1871624251L, 
726649119L, -674232154L, -757519824L, -611880749L, 451411697L, 
1167498784L, 344074430L, -774132535L, -1278343557L, -338656182L, 
-1106490196L, 2115378095L, 268296037L, -1057229540L, 944790898L, 
-770075155L, -718456281L, -780969394L, 1439446100L, -1758152046L, 
1860803360L, -1454964276L, 892021728L, 1534662338L, 618706472L, 
-1109832020L, 121288380L, -165371406L, -1376044752L, -1443821404L, 
-1136661832L, -1744621734L, -607994480L, 1003294076L, 128594196L, 
808115522L, 2080236288L, 1355475900L, 163078992L, 491738402L, 
-1815779000L, -872330228L, -1134080932L, 1335293138L, 1139659776L, 
-1447680172L, 368642536L, -2084578310L, 645319376L, -1945656020L, 
-730153740L, -1746689134L, -1312125856L, 2009323788L, -1564595488L, 
-1367485086L, -1337960792L, -213438900L, 1254011772L, 610090802L, 
1704600176L, -63893180L, -977868840L, 117552250L, -1540613584L, 
-145822020L, -938639532L, -857125310L, 2069905312L, -776078148L, 
-944098864L, -314928798L, 914270184L, 1281309644L, 396119804L, 
-1066630158L, 1125692544L, -1304362764L, 927883400L, 1285968186L, 
1917662032L, 640794732L, -2098242412L, -146250414L, 1811198688L, 
-1584173876L, 1505795232L, 1362238850L, 532024552L, 1317842156L, 
-1122473924L, 704285810L, -877527952L, -666455516L, 1126422072L, 
-61189606L, -655993648L, 294594108L, -768155500L, 1998974082L, 
53496768L, 638169404L, 1937206928L, 1981973538L, -1870745592L, 
755643532L, -429046756L, -877105006L, 1852040320L, -1410524140L, 
-1600432344L, 172410938L, 1432565968L, -920520596L, 1401360052L, 
295398418L, -1593551904L, 1855153100L, 1445667296L, -329897950L, 
-2028870360L, 49459468L, -1873632708L, -1642268430L, 1689910512L, 
-979609788L, 1187454040L, -377850310L, -740536080L, 1505228476L, 
-2002660844L, 1044872386L, 513018464L, -1733631620L, -1267774256L, 
1869637666L, 1886090408L, -857661044L, 835391932L, -261485518L, 
1363695424L, 472928948L, -1734394552L, -1616914758L, -721190128L, 
502634732L, 1410935892L, -422950126L, -1344492128L, 997832652L, 
809677664L, 1660968130L, 333660712L, -1411596884L, 182799164L, 
-348921998L, 815595312L, -1925938140L, 2045919800L, 577783514L, 
-1337457776L, -1263527684L, -1893080684L, -443014334L, -76441088L, 
-1746770884L, 1286793808L, -1975450462L, 1570633544L, -1815331316L, 
2045017820L, -290029486L, -595820416L, 1749752660L, -1492868120L, 
-2115738886L, -360702128L, -1400061908L, 1702244724L, 719164050L, 
2003493472L, -2052258676L, -1801718432L, 1100499170L, -621674328L, 
-1825117492L, 2074633724L, 1335230386L, -997566480L, 1908491716L, 
-939943592L, 1603548538L, 1855382960L, 1134586300L, 891058772L, 
430394818L, -1684063456L, 373689276L, 1953547344L, -1854332958L, 
-836295832L, 1978001228L, 1739717372L, 1044396786L, 564490496L, 
241816692L, 2146949640L, 1087138234L, -934626864L, 860335084L, 
1725925268L, 984952018L, 828622944L, 1037148492L, 1584129824L, 
-920081790L, 1650772456L, 386175980L, 647224892L, 328702578L, 
743006064L, -1315040988L, -1507537352L, -107553382L, -329681456L, 
1429803452L, 1008755348L, 833402498L, 463311296L, -1478262980L, 
698119824L, 793531170L, -54606328L, 1537521420L, 2137585564L, 
-1540218862L, 1537155968L, -1088200940L, -1921677016L, 85762106L, 
1889179600L, 590080620L, 174918708L, -1297771374L, 7601504L, 
2054151088L, 1734715353L, 1994166603L, 1847883708L, 1765781482L, 
-1990354913L, 1266089641L, 559309374L, 737095468L, 910449021L, 
-2121971769L, 2085392048L, -342926130L, 1765246459L, -2076533155L, 
-1248417398L, 2115826088L, -1761636591L, 1883368963L, -1233073916L, 
362650930L, -1291957369L, -54388271L, 1867256662L, -623421132L, 
-1746851771L, 339488911L, 1350622728L, 1356399142L, 497517267L, 
1209373749L, 88509330L, -233112736L, 1755629257L, -1602176517L, 
-66345524L, -128015270L, 1287575887L, -2014230503L, 475663918L, 
853718076L, -1476034131L, 361697879L, 564734880L, -850800386L, 
-1084550197L, -713570867L, 2125031674L, -471631048L, 501151969L, 
1864973203L, 1062518708L, -359430974L, 1234946263L, 1563063585L, 
-1052871450L, 1833607652L, 775789845L, 1308185919L, 777963736L, 
1252907254L, -1919302909L, 1637749829L, -14799774L, -1092473904L, 
324673785L, 695552555L, 739923356L, 1350987210L, -247390529L, 
1537343497L, -2112708322L, -1822022452L, -2130055011L, -1753184985L, 
290494288L, 148953582L, 126185115L, -33758979L, -852359190L, 
-906313592L, -1373733327L, 657381539L, -228979548L, 766160210L, 
-782063449L, -1139868431L, -1563571146L, -178068268L, 2110072933L, 
1212060271L, 1546964136L, -1620448762L, -714822861L, 1703362453L, 
303307762L, 1745907776L, -1608199767L, -760479717L, -913166868L, 
2027458042L, 1369881135L, -1873086919L, -774037426L, -268777316L, 
1995688461L, 1587979895L, 1638907520L, 602889438L, 2105196459L, 
-996347731L, -2013436646L, -946028072L, -2030437823L, 1755099379L, 
816549652L, -679385438L, -1104601161L, 663661441L, 997742982L, 
1259158468L, 1647486709L, 735638367L, -411099720L, -2082069354L, 
112237219L, -1314751387L, -874862206L, 838038256L, 501409433L, 
-44701045L, -171496580L, 1506807210L, -956963361L, 1778479977L, 
1966519166L, 1487088492L, 666132413L, 1236675719L, 1409212528L, 
1235782414L, -447207365L, -644765923L, 1485325770L, -855294232L, 
408211921L, -777854013L, 211860676L, 1236534130L, 449371719L, 
81559057L, 1544097686L, 176354036L, 500855941L, -776214065L, 
-1558905656L, -777622170L, 193351571L, -1097264907L, 1476761810L, 
1417946912L, -746125559L, 1167732795L, 189244812L, -129246182L, 
-2032758641L, -1600690215L, 1524353134L, 2037140476L, -1716712429L
)
