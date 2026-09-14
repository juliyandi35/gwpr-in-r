# Import libraries
library(sf)
library(tidyverse)
library(GWmodel)    # to undertake the GWR
library(tmap)       # for mapping
library(spdep)
library(ggplot2)
library(caret)
library(ggpubr)
library(gwrr)
library(readxl)
library(sp)
library(dplyr)

DATA_OK <- read_excel("DATA OK.xlsx")
head(DATA_OK)
summary(DATA_OK)

# Menghapus kata "Kabupaten " dari kolom Kab
DATA_OK$Kota <- gsub("^Kabupaten ", "", DATA_OK$Kota)
DATA_OK <- DATA_OK %>%
  arrange(Kota) %>%        # Urutkan berdasarkan Kota
  mutate(ID = row_number())  # Tambahkan kolom ID berdasarkan urutan

#Eksplorasi Data di Excel

# Import data peta
Jatim_map <- read_sf('RBI_50K_2023_Jawa Timur.x26272/RBI_50K_2023_Jawa Timur.shp')
names(Jatim_map)
Jatim_map <- Jatim_map[,c(1,26)]
colnames(Jatim_map) <- c("Kota","geometry")
Jatim_map <- Jatim_map %>%
  arrange(Kota) %>%        # Urutkan berdasarkan Kota
  mutate(ID = row_number())  # Tambahkan kolom ID berdasarkan urutan
Dataset.sdf <- as(st_zm(Jatim_map),"Spatial")

# Gabungkan data peta dan data tabular/excel
Data <- merge(Jatim_map,DATA_OK,by="Kota")
Data <- st_as_sf(Data)
Data <- st_zm(Data)

# Gambaran Data Y
ggplot(Data) + geom_sf(aes(fill = PK)) +
  scale_fill_gradient2(
    midpoint = mean(Data$PK), low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  geom_text(
    aes(label = Kota, x = coordinates(as(Data,"Spatial"))[,1], y = coordinates(as(Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = FALSE
  )+ggtitle("Sebaran Data PK di Jawa Tengah")+xlab("Longitude")+ylab("Latitude")
  

#Pembentukan Model Regresi Data Panel
library(plm)

#Model Pengaruh Acak (REM)
modelpanel1<-plm(PK~IPM+RPP+TPT+UM+GR,DATA_OK,model = "random", index = c ("Kota","Tahun"))
summary(modelpanel1)

modelpanel1_gls<-plm(PK~IPM+RPP+TPT+UM+GR,data = DATA_OK,model = "random", index = c ("Kota","Tahun"), effect = "twoways", random.method = "nerlove")
summary(modelpanel1_gls)

#Uji kebaikan REM
#Uji pengaruh individu
plmtest(modelpanel1_gls,type = "bp", effect="individu")
#Uji pengaruh waktu
plmtest(modelpanel1_gls,type = "bp", effect="time")
#Uji pengaruh twoways
plmtest(modelpanel1_gls,type = "bp", effect="twoways")


#Model Pengaruh Tetap (FEM)#
modelpanel2<-plm(PK~IPM+RPP+TPT+UM+GR,DATA_OK,model = "within", index = c ("Kota","Tahun"))
summary(modelpanel2)

#FEM INDIVIDU#
fem_individu <- plm(PK~IPM+RPP+TPT+UM+GR, index = c("Kota", "Tahun"), model = "within", effect= "individual", data = DATA_OK)
summary(fem_individu)
#Nilai diatas adalah nilai pengaruh konstan dari masing-masing individu yang pada model dapat dituliskan sebagai fi#

#FEM WAKTU#
fem_time <- plm(PK~IPM+RPP+TPT+UM+GR, index = c("Kota", "Tahun"), model = "within", effect= "time", data = DATA_OK)
summary(fem_time)
#Nilai diatas adalah nilai pengaruh konstan dari masing-masing waktu yang pada model dapat dituliskan sebagai λt#

#FEM DUA ARAH#
fem_twoways <- plm(PK~IPM+RPP+TPT+UM+GR, index = c("Kota", "Tahun"), model = "within", effect= "twoways", data = DATA_OK)
summary(fem_twoways)
data.frame(summary(fixef(fem_twoways, effect="twoways")))
#FEM DUA ARAH TIDAK LAYAK DIGUNAKAN KARENA NILAI P-VALUE >5% #

#Nilai diatas adalah nilai pengaruh konstan dari gabungan masing-masing individu dan waktu yang pada model dapat dituliskan sebagai fi+λt#


#uji Pengaruh Individu#
plmtest(fem_individu, type = "bp", effect = "individual")
#Uji Pengaruh Waktu#
plmtest(fem_time, type = "bp", effect = "time")
#Uji Pengaruh Dua Arah#
#plmtest(fem_twoways, type = "bp", effect = "twoways")  #p-value > 5%

#Karena FEM individu dan FEM time tolak H0 maka perlu dibandingkan untuk mencari model terbaik#

# Pakai syntax di bawah untuk install package nya, pastikan punya library devtools
# library(devtools)
# devtools::install_github(repo = "https://github.com/MichaelChaoLi-cpu/GWPR.light")

# Menentukan fungsi pembobot spasial terbaik
#adaptive bisquare

library(GWPR.light)
bw <- bw.GWPR(formula = PK~IPM+RPP+TPT+UM+GR, data = Dataset.df,
              index = c("Kota", "Tahun"), SDF = Dataset.sdf, adaptive = TRUE,
              effect = "twoways", model = "within",
              kernel = "bisquare", longlat = FALSE,approach = "CV")

library(GWPR.light)
result <- GWPR(bw = bw, formula = PK~IPM+RPP+TPT+UM+GR, data = Dataset.df,
               index = c("ID", "Tahun"), SDF = Dataset.sdf, adaptive = TRUE,
               effect = "twoways", model = "within",
               kernel = "bisquare", longlat = FALSE)

# Define the Exponential kernel function
library(spdep)
Dataset <- st_as_sf(st_zm(Dataset))
x = coordinates(as(Dataset,"Spatial"))[,1]
y = coordinates(as(Dataset,"Spatial"))[,2]
coords <-cbind(x,y)
jarak <-as.matrix(1/dist(coords))

bisquare_kernel <- function(dists, bw) {
  ifelse(dists < bw, (1-(dists/bw)^2)^2, 0)
}
weights <- bisquare_kernel(jarak, bw)
# writexl::write_xlsx(weights,"Pembobot Kernel Bisquare.xlsx")

# Perbandingan hasil
Compare2 <- data.frame(R.Squared = c(result$R2,summary(fem_twoways)$r.squared))
rownames(Compare2) <- c("GWPR","FEM","Adjusted R FEM")
Compare2

#Pembuatan jarak euclidean
n <- 34 #jumlah wilayah
U <- Dataset$LONG #data Longitude
V <- Dataset$LAT #data Latitude
d <- matrix(0,n,n)
for (i in 1:n) {
  for (j in 1:n) {
    d[i,j] <- sqrt(((U[i]-U[j])^2)+((U[i]-U[j])^2))
  }
}
d
# writexl::write_xlsx(data.frame(d),"Hasil_Jarak_Euclidean.xlsx")



############################TIDAK TERPAKAI#########################################
#Perbandingan kebaikan pengaruh individu dan pengaruh waktu#
#Nilai Kebaikan Model
#Sum Squared Error
dsse <- data.frame(Individu=sum(fem_individu$residuals^2),Time=sum(fem_time$residuals^2))

#AIC#
lsdv_ind <- lm(PK~IPM+RPP+TPT+UM+GR+Kota DATA_OK)
lsdv_time <- lm(PK~IPM+RPP+TPT+UM+GR + Kota + Tahun, DATA_OK)
daic <- data.frame(Individu=AIC(lsdv_ind),Time=AIC(lsdv_time))


#R-squared
ind <- summary(fem_individu)
time <- summary(fem_time)
drsq <- data.frame(Individu=ind$r.squared,Time =time$r.squared)

#Perbandingan#
compare <- t(rbind(dsse,daic,drsq))
colnames(compare) <- c("SSE","AIC","R-Squared","Adj R-Squared")
compare
############################TIDAK TERPAKAI#########################################



#Model Pengaruh Gabungan (CEM)
modelpanel3<-plm(PK~IPM+RPP+TPT+UM+GR,DATA_OK,model = "pooling", index = c ("Kota","Tahun"))
summary(modelpanel3)


#Uji Chow
pFtest(fem_individu,modelpanel3)

#Uji Chow
pooltest(modelpanel3,fem_individu)
#Model pengaruh tetap terpilih (fem pengaruh individu terpilih)


#Uji Hausman
phtest(fem_individu,modelpanel1_gls)
#Model Pengaruh acak terpilih
summary(modelpanel1_gls)
modelpanel1_gls$coefficients
#note +2 dikali 100


#Pengujian Asumsi untuk model data panel
#Uji Normalitas#
ks.test(modelpanel1_gls$residuals, "pnorm",mean=mean(modelpanel1_gls$residuals),sd=sd(modelpanel1_gls$residuals))

#Uji Autokorelasi
#pdwtest(modelpanel1)

#Multikol disimpen di atas#
#Uji Multikol#
library(car)
RegBerganda=lm(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK)
Check_Multikol=as.data.frame(t(vif(RegBerganda)))
Check_Multikol

#Uji Heterokedastisitas (Breusch-Pagan Test)
bptest(modelpanel1_gls)
#Tolak H0 artinya terdapat heterokedastisitas spasial


#Uji Kehomogenan
#library(lmtest)
#bptest(modelpanel1)


#residu_panel<- modelpanel1$residuals
#residu_panel
#residu_panel2<-as.vector(residu_panel)
#residu_panel2


#Uji Asumsi Model#
#Uji Normalitas
#kolmogorov Smirnov
#library(nortest)
#lillie.test(residu_panel2)
#Sisaan berdistribusi normal

#Uji Multikolinearitas
#Nilai VIF
#Reglin=lm(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK)
#summary(Reglin)
#library(car)
#multikol=as.data.frame(vif(Reglin))
#multikol
#Tidak adanya multikolinearitas


#Uji Autokorelasi
#durbinWatsonTest(residu_panel2)
#nilai du= 1.8063, dl=1,6776, nilai d = 1.559211


#Uji Heterokedastisitas (Breusch-Pagan Test)
#library(lmtest)
#yresidu =abs(residu_panel2)
#bptest(modelpanel1,studentize = FALSE)
#Tolak H0 artinya terdapat heterokedastisitas spasial


#Model GWPR#
library(spgwr)
library(GWmodel)
head(DATA_OK,38)
getOption("max.print")
library(sp)

#Input Latitude Longitude
coordinates(DATA_OK)=4:5
class(DATA_OK)
head(DATA_OK,38)


#MENENTUKAN Fungsi pembobot spasial terbaik#

#Menghitung jarak euclidean
coordinates(DATA_OK)
coords <- coordinates(DATA_OK)
View(coords)

jarak_euclidean <- dist(coords,method = "euclidean")
jarak_matrix<-as.matrix(jarak_euclidean)
View(jarak_matrix)
write.csv(jarak_matrix, "jarak_euclidean.csv")


#Menentukan Fungsi Pembobot Spasial Terbaik

#Adaptive Bisquare
#result<-bw.gwr(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,approach = "CV")
bwd.GWPR.bisquare.ad <- bw.gwr(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,approach = "CV",kernel = "bisquare",adaptive = T)
hasil.GWPR.bisquare.ad<-gwr.basic(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,bw=bwd.GWPR.bisquare.ad,kernel = "bisquare",adaptive = T)
hasil.GWPR.bisquare.ad

#Adaptive Gaussian
bwd.GWPR.gaussian.ad <- bw.gwr(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,approach = "CV",kernel = "gaussian",adaptive = T)
hasil.GWPR.gaussian.ad<-gwr.basic(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,bw=bwd.GWPR.gaussian.ad,kernel = "gaussian",adaptive = T)
hasil.GWPR.gaussian.ad

#Fixed Bisquare
bwd.GWPR.bisquare.fix <- bw.gwr(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,approach = "CV",kernel = "bisquare",adaptive = F)
hasil.GWPR.bisquare.fix<-gwr.basic(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,bw=bwd.GWPR.bisquare.fix,kernel = "bisquare",adaptive = F)
hasil.GWPR.bisquare.fix

#Fixed Gaussian
bwd.GWPR.gaussian.fix <- bw.gwr(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,approach = "CV",kernel = "gaussian",adaptive = F)
hasil.GWPR.gaussian.fix<-gwr.basic(PK~IPM+RPP+TPT+UM+GR,data=DATA_OK,bw=bwd.GWPR.gaussian.fix,kernel = "gaussian",adaptive = F)
hasil.GWPR.gaussian.fix

#Perbandingan Kernel Terbaik#
AIC<-c(hasil.GWPR.bisquare.ad$GW.diagnostic$AIC,hasil.GWPR.gaussian.ad$GW.diagnostic$AIC,hasil.GWPR.bisquare.fix$GW.diagnostic$AIC,hasil.GWPR.gaussian.fix$GW.diagnostic$AIC)
R2<- c(hasil.GWPR.bisquare.ad$GW.diagnostic$gw.R2,hasil.GWPR.gaussian.ad$GW.diagnostic$gw.R2,hasil.GWPR.bisquare.fix$GW.diagnostic$gw.R2,hasil.GWPR.gaussian.fix$GW.diagnostic$gw.R2)
AICc<-c(hasil.GWPR.bisquare.ad$GW.diagnostic$AICc,hasil.GWPR.gaussian.ad$GW.diagnostic$AICc,hasil.GWPR.bisquare.fix$GW.diagnostic$AICc,hasil.GWPR.gaussian.fix$GW.diagnostic$AICc)
R2adj<-c(hasil.GWPR.bisquare.ad$GW.diagnostic$gwR2.adj,hasil.GWPR.gaussian.ad$GW.diagnostic$gwR2.adj,hasil.GWPR.bisquare.fix$GW.diagnostic$gwR2.adj,hasil.GWPR.gaussian.fix$GW.diagnostic$gwR2.adj)
RSS <- c(hasil.GWPR.bisquare.ad$GW.diagnostic$RSS.gw,hasil.GWPR.gaussian.ad$GW.diagnostic$RSS.gw,hasil.GWPR.bisquare.fix$GW.diagnostic$RSS.gw,hasil.GWPR.gaussian.fix$GW.diagnostic$RSS.gw)


# Misalkan df memiliki kolom lat dan lon untuk koordinat geografis
coordsgwpr <- df[,c("Lat", "Long")]  # Koordinat geografis


# Model GWPR
modelgwpr <- gwr(y ~ x1 + x2, data = df, coords = coords, bandwidth = bandwidth)  # Sesuaikan bandwidth




#perbandingan kernel terbaik
tabel <- cbind(RSS,R2,R2adj,AIC,AICc)
rownames(tabel) <- c("Adaptive Bisquare", "Adaptive Gaussian", "Fixed Bisquare", "Fixed Gaussian")
data.frame(tabel)
bwd.GWPR.bisquare.ad


#class(bwd.GWPR.bisquare.ad)
#str(bwd.GWPR.bisquare.ad)
#is.atomic(bwd.GWPR.bisquare.ad)
library(sp)
library(rgdax)

#Penduga Parameter GWPR BISQUARE

library(dplyr)
parameter.GWPR <- as.data.frame(hasil.GWPR.bisquare.ad$SDF[1:38,1:7])[-c(7:7)]
parameter.GWPR%>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())

#Model GWPR
Modelgwpr <- as.data.frame(hasil.GWPR.bisquare.ad$SDF[1:38,c(1:6,5)])[-7:-9]
Modelgwpr

#Model GWPR dengan jumlah penduduk miskin per luas wilayah tingi
parameter_R2 <- as.data.frame(hasil.GWPR.bisquare.ad$SDF[1:38,c(1:6,5)])[-7:-9]
parameter_R2 %>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>%
  filter(Kota %in% c("Kabupaten Malang", "Kabupaten Jember", "Kabupaten Sampang")) %>% 
  select(Kota, everything())


#signifikasiparameter
#Uji Simultan
#Uji Kesesuaian Model
#Uji F
rss_rem <- sum(modelpanel1_gls$residuals^2)
rss_gwpr <- hasil.GWPR.bisquare.ad$GW.diagnostic$RSS.gw
df_rem <- modelpanel1_gls$df.residual
df_gwpr <- hasil.GWPR.bisquare.ad$GW.diagnostic$edf
Fhit=(rss_gwpr/df_gwpr)/(rss_rem/df_rem)
Fhit
pvalue <- pf(Fhit,df_rem,df_gwpr,lower.tail = FALSE)
pvalue

p_value<-gwr.t.adjust(hasil.GWPR.bisquare.ad)$results$p
p_value<-as.data.frame(ifelse(p_value<=0.05,"signifikan","Tidak"))[1:38,]
p_value %>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())


#R2 setiap wilayah


#######################################################################################
#Penduga Parameter GWPR ADAPTIVE BISQUARE

#Uji Simultan
#Uji Kesesuaian Model
#Uji F
rss_rem <- sum(modelpanel1$residuals^2)
rss_gwpr <- hasil.GWPR.bisquare.ad$GW.diagnostic$RSS.gw
df_rem <- modelpanel1$df.residual
df_gwpr <- hasil.GWPR.bisquare.ad$GW.diagnostic$edf
Fhit=(rss_gwpr/df_gwpr)/(rss_rem/df_rem)
Fhit
pvalue <- pf(Fhit,df_rem,df_gwpr,lower.tail = FALSE)
pvalue

#################################################
#Simpan SSR dari random effect model
#SSR_REM <-Sum(residuals(modelpanel1)^2)
#df_REM <- df.residual(modelpanel1)
#SSR dari model GWPR
#SSR_gwpr <-sum(parameter.GWPR$SDF$residual^2)
#df_gwpr <- nrow(DATA_OK)-length(parameter.GWPR$SDF@data)
#F-hitung
#df1<- df_REM - df_gwpr
#df1
#df2<-df_gwpr
#Fhit<-((SSR_REM - SSR_gwpr)/df1)/(SSR_gwpr/df2)
#alpha<-0.05
#F_crit<-qf(1-alpha, df1, df2 )
#hasil
#cat("Fhitung =",Fhit,"\n")
#cat("F-kritis =",F_crit,"\n")
#if(Fhit > F_crit){
#cat("Kesimpulan:Tolak H0 artinya model GWPR signifikan secara simultan.\n")
#} else {
#cat("Kesimpulan:Terima H0 artinya model GWPR tidak signifikan secara simultan.\n")
#}


#Uji Parsial
library(dplyr)
parameter.GWPRadapbis <- as.data.frame(hasil.GWPR.gaussian.ad$SDF[1:38,1:7])[-c(7:7)]
parameter.GWPRadapbis%>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())
summary(parameter.GWPRadapbis)

library(plm)
modelgwpr_kabblitar<-plm(PK~IPM+RPP+TPT+UM+GR,DATA_OK,model = "parameter.GWPRadapbis", index = c ("Kota","Tahun"))
summary(modelgwpr_kabblitar)
# Model GWPR
modelgwpr <- gwr(PK~IPM+RPP+TPT+UM+GR, data = df, coords = jarak_euclidean, bandwidth = )




library(GWmodel)
gwpr_result <- gwpr(formula = y ~ x1 + x2, data = panel_data, coords = coords, ...)
summary(gwpr_result$SDF)


#######################################################################################
#Penduga Parameter GWPR GAUSSIAN
library(dplyr)
parameter.GWPRgaussian <- as.data.frame(hasil.GWPR.gaussian.ad$SDF[1:38,1:7])[-c(7:7)]
parameter.GWPRgaussian%>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())

#Penduga Parameter GWPR FIX BISQUARE
library(dplyr)
parameter.GWPRfixbis <- as.data.frame(hasil.GWPR.bisquare.fix$SDF[1:38,1:7])[-c(7:7)]
parameter.GWPRfixbis%>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())

#Penduga Parameter GWPR FIX GAUSSIAN
library(dplyr)
parameter.GWPRfixgauss <- as.data.frame(hasil.GWPR.gaussian.fix$SDF[1:38,1:7])[-c(7:7)]
parameter.GWPRfixgauss%>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())


#Model GWPR dengan Kemiskinan per IPM per Luas Wilayah
#parameter.GWPR1 <- as.data.frame(hasil.GWPR.bisquare.ad$SDF[1:38,1:11])[-c(3:5)]
#parameter.GWPR1%>% 
#  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
#  filter(Kota %in% c("Kabupaten Pasuruan","Kabupaten Sumenep","Kabupaten Bangkalan")) %>%
#  select(Kota, everything())


#signifikasiparameter(pvalue)
p_value<-gwr.t.adjust(hasil.GWPR.bisquare.ad)$results$p
p_value<-as.data.frame(ifelse(p_value<=0.05,"signifikan","Tidak"))[1:38,]
p_value %>% 
  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
  select(Kota, everything())

#FILTERKOTA YANG SIGNIFIKAN PADA PEUBAH KEMISKINAN
#p_value %>% 
#  mutate(Kota=DATA_OK$Kota[1:38]) %>% 
#  select(Kota, everything())
#  filter("Signifikan" == p_value$IPM_p)
