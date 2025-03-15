# Projekt - Analiza danych jakościowych Łukasz Obrzut ---------------------
#install.packages("caret")
#install.packages("ROCR")
library(caret)
library(ROCR)
set.seed(123)

# Ładowanie danych oraz zmiana nazw kolumn ------------------------------
#setwd("")
dane <- read.csv("ThoraricSurgery.csv", header=TRUE)

# usunięcie zmiennej id(niepotrzebna w analizie)
dane <- dane[, -1]


# najpierw zamienię kolumny na inne nazwy z powodu zrozumienia
# PRE4 <-  FVC 
# PRE5 <- FEV1  
# PRE6 <- Zubrod
# PRE7 <-  Bol przed operacją 
# PRE8 <- Krwioplucie przed operacją 
# PRE9 <- 	Duszność przed operacją 
# PRE10 <- 	Kaszel przed operacją 
# PRE11 <- 	Osłabienie przed operacją 
# PRE14 <- WielGuza (4 poziomy)
# PRE19 <- Cukrzyca2 
# PRE25 <- ZawałMS
# PRE25 <- Choroba tętnic obwodowych (PAD) 
# PRE30 <- Palenie (P,K)
# PRE32 <- Astma 

# zmiana nazw kolumn
colnames(dane) <- c("DGN", "FVC", "FEV1", "Zubrod", "Bol", "Krwioplucie","Duszność", "Kaszel",
                    "Oslabienie","WielGuza", "Cukrzyca2", "ZawałMS", "PAD", "Palenie", "Astma", "Wiek", "Ryzyko1R")

# Zmiana typów zmiennych
dane$DGN <- factor(dane$DGN)
dane[,-c(2,3,16)] <- lapply(dane[,-c(2,3,16)], function(x) if(is.logical(x)) factor(x) else x)
dane$WielGuza = factor(dane$WielGuza, ordered = TRUE)
dane$Zubrod <- factor(dane$Zubrod, ordered = TRUE)

str(dane)
# Przegląd poprawności danych ---------------------------------------------
barplot(table(dane$Ryzyko1R), main = "Zgon w ciągu roku od zabiegu(TRUE-pacjent zmarł, FALSE-pacjent przeżył) ") #dane niezbalansowane
summary(dane)

### Widzimy w zmiennej FEV1 bardzo podejrzane informację. Wartość największa dla tej zmiennej jest zbyt duża,
# być może wynika to z błedu interpretacji pomiaru(można by to traktować jako procent optymalnej wartości wydychanej FEV1, jednak nie mamy takiej informacji)
# spójrzmy na histogram

hist(dane$FEV1, xlab = "FEV1", main = "histogram dla zmiennej FEV1") # widać bardzo duże wartości zmiennej FEV1 będące outlierami(taka wartość FEV1 jest nieprawdopodobna)
dane$FEV1[which(dane$FEV1>6)] #(5-6 litrów to stanadrdowa pojemność płuc dorosłego mężczyzny, widzimy ze powyżej
# 6 litrów najmniejszą wartością jest 8.56, która również jest nieprawdopodobna.

### Usuwamy te 15 obserwacji
indeksy_do_usun <- which(dane$FEV1>6)
dane = dane[-indeksy_do_usun,]
row.names(dane) <- NULL

summary(dane)

# Miary zależności -----------------------------------------------------

### Dla zmiennych kategorycznych ----------------------------------------

#### Funkcja tworząca tabelę kontyngencji i zwracająca p wartość testu Fishera  --------

pwartosc_fi <- function(x, y){
  tabela_kontyngencji= table(x, y)
  p_wartosc = fisher.test(tabela_kontyngencji, simulate=T, B=20000)$p.value
  return(p_wartosc)
}

#### Liczenie p wartości testu Fishera dla każdej pary zmiennych jakościowych  --------
zmienne_ciagle <- c("FVC", "FEV1", "Wiek")
zmienne_jakosciowe <- colnames(dane[,-c(2,3,16)])

p_macierz_fi <- matrix(0, nrow = 14, ncol = 14)
dimnames(p_macierz_fi) <- list(zmienne_jakosciowe, zmienne_jakosciowe)

for(i in zmienne_jakosciowe){
  for(j in zmienne_jakosciowe){
    p_macierz_fi[i, j] = pwartosc_fi(dane[,i], dane[,j])
  }
}

#p_macierz_fi

### Współczynniki Goodmana Kruskala ----------------------------------

# funkcja z wykładu
GK.gamma<-function(x, pr=0.95)
{
  # x is a matrix of counts.  You can use output of crosstabs or xtabs in R.
  # A matrix of counts can be formed from a data frame by using design.table.
  
  # Confidence interval calculation and output from Greg Rodd
  
  # Check for using S-PLUS and output is from crosstabs (needs >= S-PLUS 6.0)
  if(is.null(version$language) && inherits(x, "crosstabs")) { oldClass(x)<-NULL; attr(x, "marginals")<-NULL}
  
  n <- nrow(x)
  m <- ncol(x)
  pi.c<-pi.d<-matrix(0,nr=n,nc=m)
  
  row.x<-row(x)
  col.x<-col(x)
  
  for(i in 1:(n)){
    for(j in 1:(m)){
      pi.c[i, j]<-sum(x[row.x<i & col.x<j]) + sum(x[row.x>i & col.x>j])
      pi.d[i, j]<-sum(x[row.x<i & col.x>j]) + sum(x[row.x>i & col.x<j])
    }
  }
  
  C <- sum(pi.c*x)/2
  D <- sum(pi.d*x)/2
  
  psi<-2*(D*pi.c-C*pi.d)/(C+D)^2
  sigma2<-sum(x*psi^2)-sum(x*psi)^2
  
  gamma <- (C - D)/(C + D)
  pr2 <- 1 - (1 - pr)/2
  CIa <- qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + gamma
  
  list(gamma = gamma, C = C, D = D, sigma = sqrt(sigma2), Level = paste(
    100 * pr, "%", sep = ""), CI = paste(c("[", max(CIa[1], -1), 
                                           ", ", min(CIa[2], 1), "]"), collapse = ""))     
}


#### Funkcja tworząca tabelę kontyngencji i zwracająca wartość gamma  --------
gamma_GK <- function(x, y){
  tabela_kontyngencji= table(x, y)
  GK_gamma = GK.gamma(tabela_kontyngencji)$gamma
  return(GK_gamma)
}

gamma_GK_CI <- function(x, y){
  tabela_kontyngencji= table(x, y)
  GK_gamma = GK.gamma(tabela_kontyngencji)$CI
  return(GK_gamma)
}

gamma_macierz_CI <- matrix(0, nrow = 14, ncol = 14)
dimnames(gamma_macierz_CI) <- list(zmienne_jakosciowe, zmienne_jakosciowe)

#### Liczenie wartości gamma dla każdej pary zmiennych jakościowych --------
gamma_macierz_GK <- matrix(0, nrow = 14, ncol = 14)
dimnames(gamma_macierz_GK) <- list(zmienne_jakosciowe, zmienne_jakosciowe)

for(i in zmienne_jakosciowe){
  for(j in zmienne_jakosciowe){
    gamma_macierz_GK[i, j] = gamma_GK(dane[,i], dane[,j])
    gamma_macierz_CI[i, j] = gamma_GK_CI(dane[,i], dane[,j])
  }
}

#gamma_macierz_GK



# TAU - funkcja z wykładu 
GK.tau <- function(dat)
{ N <- sum(dat);

dat.rows <- nrow(dat);
dat.cols <- ncol(dat);
max.col <- sum.col <- L.col <- matrix(,dat.cols);
max.row <- sum.row <- L.row <- matrix(,dat.rows);
for(i in 1:dat.cols)
{ sum.col[i] <- sum(dat[,i]); max.col[i] <- max(dat[,i]); }
for(i in 1:dat.rows)
{ sum.row[i] <- sum(dat[i,]); max.row[i] <- max(dat[i,]); }

max.row.margin <- max(apply(dat,1,sum));   max.col.margin <- max(apply(dat,2,sum));

# Goodman-Kruskal tau (raws=indep.vars, cols=dep.vars)
n.err.unconditional <- N^2;
for(i in 1:dat.rows)
  n.err.unconditional <- n.err.unconditional-N*sum(dat[i,]^2/sum.row[i]);   
n.err.conditional <- N^2-sum(sum.col^2);   
tau <- 1-(n.err.unconditional/n.err.conditional);

v <- n.err.unconditional/(N^2);
d <- n.err.conditional/(N^2);
f <- d*(v+1)-2*v;

var.tau.CR <- 0;
for(i in 1:dat.rows)
  for(j in 1:dat.cols)
    var.tau.CR <- var.tau.CR + dat[i,j]*(-2*v*(sum.col[j]/N)+d*((2*dat[i,j]/sum.row[i])-sum((dat[i,]/sum.row[i])^2))-f)^2/(N^2*d^4);
ASE <- sqrt(var.tau.CR);

U.tau.CR <- (N-1)*(dat.cols-1)*tau; 
# Chi-squared approximation for H0 according to Margolin & Light JASA 1974, 755-764, 
# see also Liebetrau 1983   
p.value <- pchisq(U.tau.CR,df=(dat.rows-1)*(dat.cols-1),lower=FALSE); 

data.frame(tau, p.value, ASE);  
}

#### Funkcja tworząca tabelę kontyngencji i zwracająca p wartość dla tau  --------
tau_GK <- function(x, y){
  tabela_kontyngencji= table(x, y)
  GK_tau = GK.tau(tabela_kontyngencji)$p.value
  return(GK_tau)
}


#### Liczenie p wartości testu dla każdej pary zmiennych jakościowych  --------
tau_macierz_GK <- matrix(0, nrow = 14, ncol = 14)
dimnames(tau_macierz_GK) <- list(zmienne_jakosciowe, zmienne_jakosciowe)

for(i in zmienne_jakosciowe){
  for(j in zmienne_jakosciowe){
    tau_macierz_GK[i, j] = tau_GK(dane[,i], dane[,j])
  }
}



## Usuwam Astma oraz ZawałMS ze względu na tylko dwie wartości True
dane <- dane[,-c(12,15)]
colnames(dane)
summary(dane)
#### podział na zbiór testowy i treningowy w stosunku 85% treningowy, 15% testowy -------------------------------
# używam funkcji createDataPartition, ponieważ chcemy aby nasze klasy były zrównoważone zarówno w zbiorze testowym i treningowym

indeksy <- createDataPartition(dane$Ryzyko1R, p = 0.85, list = FALSE)

dane_treningowe <- dane[indeksy, ]  
dane_testowe <- dane[-indeksy, ] 

summary(dane_treningowe$Ryzyko1R)


####Backward -------------------------------------------------------------

## Model 1 -----------------------------------------------------------------

model_pelny <- glm(Ryzyko1R ~ ., family = binomial, data = dane_treningowe)
summary(model_pelny) # dużo zmiennych nieistotnych

backward <- step(model_pelny, list(lower = ~ 1, upper =formula(model_pelny)), trace = F, direction = "backward")
backward$anova
coef(backward)
names(coef(backward))[-1]

model_1 <- glm(Ryzyko1R ~ DGN + FEV1 + Bol + Duszność + WielGuza + Palenie, family = binomial, data = dane_treningowe)
summary(model_1) #AIC = 316.99, problem ze zmienna DGN

##### Forward -------------------------------------------------------------
model_prosty = glm(Ryzyko1R ~ 1, family = binomial, data = dane_treningowe)
forward <- step(model_prosty, 
                scope = list(lower = ~ 1, upper = model_pelny), 
                direction = "forward")
forward$anova
coef(forward)
names(coef(forward))[-1]
names(coef(backward))[-1]

#Oba modele wskazały te same zmienne w kontekscie minimalizacji AIC

##### Sprawdzenie jak zachowuje się Gmean w zależności od ucięcia -------------

klasa_predykcje <- prediction(predict(model_1, type = "response"), dane_treningowe$Ryzyko1R)

# TPR i TNR
perf_tpr <- performance(klasa_predykcje, "tpr")
perf_tnr <- performance(klasa_predykcje, "tnr")

TPR <- perf_tpr@y.values[[1]]
TNR <- perf_tnr@y.values[[1]]

G_mean_x <- sqrt(TPR * TNR)

plot(perf_tpr@x.values[[1]],G_mean_x ,
     xlab = "Próg", ylab = "Gmean",
     main = "Wykres Gmean w zależności od progu")
##### Próg maksymalizujący gmean-----------------

indeks_maksymalnego <- which.max(G_mean_x)
(max_próg <- perf_tnr@x.values[[1]][indeks_maksymalnego])

##### Wyniki Gmean i macierz błedu dla tego poziomu ----------------------
przewidywane_prawd <- predict(object=model_1,type="response", newdat = dane_treningowe)

przewidziana_klasa <- ifelse(przewidywane_prawd > 0.141677, 1, 0)

macierz_bledu <- table(przewidziana_klasa, dane_treningowe$Ryzyko1R)

rownames(macierz_bledu) <- c("przewidziane przezycie",
                             "przewidziana smierc")
colnames(macierz_bledu) <- c("obserwowane_przezycie",
                             "obserwowana_smierc")
print(macierz_bledu)  

tnr = macierz_bledu[1,1]/(macierz_bledu[1,1]+macierz_bledu[2,1])
tnr

tpr = macierz_bledu[2,2]/(macierz_bledu[1,2]+macierz_bledu[2,2])
tpr

G_mean = sqrt(tpr*tnr)
G_mean # 0.719, jednak model ten nie wydaje się najlepszy

res_pearson_m1 <- residuals(model_1, type = "pearson" )
res_m1 <-  residuals(model_1)

plot(res_pearson_m1, main = "Reszty Pearsonowskie")
plot(res_m1)

### Model 2--------------------------------------------

### Spróbuję dodać zmienną z interakcją Palenie:WielGuza oraz WielGuza:Krwioplucie

model_pelny_m2 <- glm(Ryzyko1R ~ . + WielGuza:Palenie + WielGuza:Krwioplucie, family = binomial, data = dane_treningowe)
summary(model_pelny_m2) # dużo zmiennych nieistotnych

forward <- step(model_prosty, 
                scope = list(lower = ~ 1, upper = model_pelny_m2), 
                direction = "forward")
forward$anova
coef(forward)
names(coef(forward))[-1]

## dodatkowo wskazuje wielGuza:Palenie
model_2 <- glm(Ryzyko1R ~ WielGuza + DGN + FEV1 + Bol + Palenie + Duszność + WielGuza:Palenie, family = binomial, data = dane_treningowe)
summary(model_2) # AIC - 307.3872, zmniejsza się AIC, ale model wydaje się być nienajlepszy
przewidywane_prawd <- predict(object=model_2,type="response", newdat = dane_treningowe)

##### Sprawdzenie jak zachowuje się Gmean w zależności od ucięcia -------------

klasa_predykcje <- prediction(predict(model_2, type = "response"), dane_treningowe$Ryzyko1R)

# TPR i TNR
perf_tpr <- performance(klasa_predykcje, "tpr")
perf_tnr <- performance(klasa_predykcje, "tnr")

TPR <- perf_tpr@y.values[[1]]
TNR <- perf_tnr@y.values[[1]]

G_mean_x <- sqrt(TPR * TNR)

plot(perf_tpr@x.values[[1]],G_mean_x ,
     xlab = "Próg", ylab = "Gmean",
     main = "Wykres Gmean w zależności od progu")

##### Próg maksymalizujący gmean-----------------

indeks_maksymalnego <- which.max(G_mean_x)
(max_próg <- perf_tnr@x.values[[1]][indeks_maksymalnego])

res_pearson_m2 <- residuals(model_2, type = "pearson" )
res_m2 <-  residuals(model_2)

plot(res_pearson_m2, main = "Reszty Pearsonowskie")
plot(res_m2)

##### Wyniki Gmean i macierz błedu dla tego poziomu ----------------------


przewidziana_klasa <- ifelse(przewidywane_prawd > 0.165, 1, 0)

macierz_bledu <- table(przewidziana_klasa, dane_treningowe$Ryzyko1R)

rownames(macierz_bledu) <- c("przewidziane przezycie",
                             "przewidziana smierc")
colnames(macierz_bledu) <- c("obserwowane_przezycie",
                             "obserwowana_smierc")
print(macierz_bledu) 

tnr = macierz_bledu[1,1]/(macierz_bledu[1,1]+macierz_bledu[2,1])
tnr

tpr = macierz_bledu[2,2]/(macierz_bledu[1,2]+macierz_bledu[2,2])
tpr

G_mean = sqrt(tpr*tnr)
G_mean # 0.715

res_pearson_m1 <- residuals(model_1, type = "pearson" )
res_m1 <-  residuals(model_1)

plot(res_pearson_m1, main = "Reszty Pearsonowskie")
plot(res_m1)

### Model 3--------------------------------------------

### Spróbuję usunąc zmienną DGN z rozważań modelu pełnego, gdyż daje nam w modelach wysokie p wartości(teoretycznie nieistotna),
# jednak w kontekscie minimalizacji AIC dobrze się sprawdza
# Usunę zmienną DGN z rozważanego modelu

model_pelny_m3 <- glm(Ryzyko1R ~ . - DGN, family = binomial, data = dane_treningowe)
summary(model_pelny_m3) # dużo zmiennych nieistotnych

forward <- step(model_prosty, 
                scope = list(lower = ~ 1, upper = model_pelny_m3), 
                direction = "forward")
forward$anova
coef(forward)
names(coef(forward))[-1]

## model3
model_3 <- glm(Ryzyko1R ~ WielGuza + FEV1 + Bol + Palenie + Duszność, family = binomial, data = dane_treningowe)
summary(model_3) # AIC - 324, model ma większe AIC, jednka zmienne w modelu wydają się najbardziej sensowne
przewidywane_prawd <- predict(object=model_3, type="response", newdat = dane_treningowe)

library(regclass)
#VIF(model_3)
#VIF(model_1)
#VIF(model_2)
##### Sprawdzenie jak zachowuje się Gmean w zależności od ucięcia -------------

klasa_predykcje <- prediction(predict(model_3, type = "response"), dane_treningowe$Ryzyko1R)

# TPR i TNR
perf_tpr <- performance(klasa_predykcje, "tpr")
perf_tnr <- performance(klasa_predykcje, "tnr")

TPR <- perf_tpr@y.values[[1]]
TNR <- perf_tnr@y.values[[1]]

G_mean_x <- sqrt(TPR * TNR)

plot(perf_tpr@x.values[[1]],G_mean_x ,
     xlab = "Próg", ylab = "Gmean",
     main = "Wykres Gmean w zależności od progu")

##### Próg maksymalizujący gmean-----------------

indeks_maksymalnego <- which.max(G_mean_x)
(max_próg <- perf_tnr@x.values[[1]][indeks_maksymalnego])


res_pearson_m3 <- residuals(model_3, type = "pearson" )
res_m3 <-  residuals(model_3)

plot(res_pearson_m3, main = "Reszty Pearsonowskie")
plot(res_m3)
##### Wyniki Gmean i macierz błedu dla tego poziomu ----------------------

przewidziana_klasa <- ifelse(przewidywane_prawd > 0.16, 1, 0)

macierz_bledu <- table(przewidziana_klasa, dane_treningowe$Ryzyko1R)

rownames(macierz_bledu) <- c("przewidziane przezycie",
                             "przewidziana smierc")
colnames(macierz_bledu) <- c("obserwowane_przezycie",
                             "obserwowana_smierc")
print(macierz_bledu) 

tnr = macierz_bledu[1,1]/(macierz_bledu[1,1]+macierz_bledu[2,1])
tnr

tpr = macierz_bledu[2,2]/(macierz_bledu[1,2]+macierz_bledu[2,2])
tpr

G_mean = sqrt(tpr*tnr)
G_mean # 0.68

### Jako, że wskażnik gmean pokazuje podobny trend dla tych samych progów decyduje się na sprawdzenie
### za pomocą testu ilorazu wiarygodności czy model 2 jest istotnie lepszy(wzgledem kryterium dewiancji) niz 1 oraz czy model 1 jest istotnie lepszy niż 2 
anova(model_1, model_2, test = "Chisq")
anova(model_3, model_1, test = "Chisq")
anova(model_3, model_2, test = "Chisq")
anova(model_3, model_1, model_2, test = "Chisq")

##Sprawdzenie na zbiorze testowym------------------------------------------------------------

# wybieram model nr 2
przewidywane_prawd_test <- predict(object=model_2, type="response", newdata = dane_testowe)

przewidziana_klasa_test <- ifelse(przewidywane_prawd_test > 0.16, 1, 0) # wybieram 0.15 ponieważ z wykresu Gmean widzimy, 
# że do 0.16 ma tendencje rosącą, jednak przed 0.16 widzimy znaczną poprawe

macierz_bledu_test <- table(przewidziana_klasa_test, dane_testowe$Ryzyko1R)

rownames(macierz_bledu_test) <- c("przewidziane przezycie",
                                  "przewidziana smierc")
colnames(macierz_bledu_test) <- c("obserwowane_przezycie",
                                  "obserwowana_smierc")
print(macierz_bledu_test) 

tnr_test = macierz_bledu_test[1,1]/(macierz_bledu_test[1,1]+macierz_bledu_test[2,1])
tnr_test

tpr_test = macierz_bledu_test[2,2]/(macierz_bledu_test[1,2]+macierz_bledu_test[2,2])
tpr_test

G_mean_test = sqrt(tpr_test*tnr_test)
G_mean_test # 0.592