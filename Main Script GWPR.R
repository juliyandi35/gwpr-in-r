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
library(car)
library(broom)
library(lmtest)
# library(devtools)
# devtools::install_github(repo = "https://github.com/MichaelChaoLi-cpu/GWPR.light")
library(GWPR.light)

Dataset <- read_excel("DATA OK.xlsx")
names(Dataset)
colnames(Dataset) <- c("ID","Kota","Tahun","Lat","Long","PK","IPM","RPP","TPT","UM","GR")
head(Dataset)
summary(Dataset)

# Menghapus kata "Kabupaten " dari kolom Kab
Dataset$Kota <- gsub("^Kabupaten ", "", Dataset$Kota)

#Eksplorasi Data di Excel

# Import data peta
Jatim_map <- read_sf('RBI_50K_2023_Jawa Timur.x26272/RBI_50K_2023_Jawa Timur.shp')
names(Jatim_map)
Jatim_map <- Jatim_map[,c(1,26)]
colnames(Jatim_map) <- c("Kota","geometry")
Jatim_map<-merge(Jatim_map,Dataset[,c(1,2)],by="Kota")
Dataset.sdf <- as(st_zm(Jatim_map),"Spatial")

# Gabungkan data peta dan data tabular/excel
Data <- merge(Jatim_map,Dataset,by="ID")
Data <- st_as_sf(Data)
Data <- st_zm(Data)
Data <- Data[,-2]
colnames(Data) <- c("ID","Kota","Tahun","Lat","Long","PK","IPM","RPP","TPT","UM","GR","geometry")

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

# Moran test
# Perbaiki hanya yang tidak valid
Data_Moran <- Data
Data_Moran$geometry[!st_is_valid(Data_Moran)] <- st_make_valid(Data_Moran$geometry[!st_is_valid(Data_Moran)])

# Tangani geometri kompleks yang pecah jadi multipolygon
Data_Moran <- Data_Moran %>% 
  st_cast("MULTIPOLYGON") %>%  # ubah semua ke multipolygon
  group_by(ID) %>%             # gabungkan lagi berdasarkan ID asli
  summarise(across(everything(), first), .groups = "drop")  # ambil info non-geometry

# Sekarang jumlah baris kembali seperti semula
nrow(Data_Moran)  # Harusnya sama seperti sebelum validasi

# Buat objek tetangga & spatial weights
sf::sf_use_s2(FALSE)
nb <- poly2nb(Data_Moran)
lw <- nb2listw(nb, style = "W")

# Uji Moran
moran.test(Data_Moran$PK, lw)

# Uji Multikolinieritas
check_model <- lm(PK~IPM+RPP+TPT+UM+GR,data=Dataset)
data.frame(t(vif(check_model)))

# Model Common Effect
library(plm)
cem <- plm(PK~IPM+RPP+TPT+UM+GR,data=Dataset, model = "pooling")
summary(cem)

# Model Fixed Effect
fem <- plm(PK~IPM+RPP+TPT+UM+GR,data=Dataset, index = c("Kota", "Tahun"), model = "within", effect= "individual")
summary(fem)
summary(fixef(fem, effect="individual"))

# Model FEM dengan waktu
fem_time <- plm(PK~IPM+RPP+TPT+UM+GR,data=Dataset, index = c("Kota", "Tahun"), model = "within", effect= "time")
summary(fem_time)
summary(fixef(fem_time, effect="time"))

# Model FEM 2 arah
fem_twoways <- plm(PK~IPM+RPP+TPT+UM+GR,data=Dataset, index = c("Kota", "Tahun"), model = "within", effect= "twoways")
summary(fem_twoways)
data.frame(summary(fixef(fem_twoways, effect="twoways")))

# Uji Signifikasi Pengaruh Individu / Waktu / Two ways
# Uji pengaruh individu
plmtest(fem_twoways, type = "bp", effect = "individual")

# Uji pengaruh waktu
plmtest(fem_twoways, type = "bp", effect = "time")

# Uji pengaruh twoways
plmtest(fem_twoways, type = "bp", effect = "twoways")

# Nilai Kebaikan Model
# Sum Squared Error
dsse <- data.frame(Individu=sum(fem$residuals^2),Time=sum(fem_time$residuals^2),Twoways=sum(fem_twoways$residuals^2))

# AIC
lsdv_ind <- lm(PK~IPM+RPP+TPT+UM+GR + Kota,data=Dataset)
lsdv_time <- lm(PK~IPM+RPP+TPT+UM+GR + Tahun,data=Dataset)
lsdv_twoways <- lm(PK~IPM+RPP+TPT+UM+GR + Kota + Tahun,data=Dataset)

daic <- data.frame(Individu=AIC(lsdv_ind),Time=AIC(lsdv_time),Twoways=AIC(lsdv_twoways))

# MAPE
mape <- function(actual, forecast) {
  mean(abs((actual - forecast) / actual)) * 100
}
dmape <- data.frame(Individu=mape(Dataset$PK,predict(fem)),
                    Time=mape(Dataset$PK,predict(fem_time)),
                    Twoways=mape(Dataset$PK,predict(fem_twoways)))
# BIC
dbic <- data.frame(Individu=BIC(lsdv_ind),Time=BIC(lsdv_time),Twoways=BIC(lsdv_twoways))

# R-squared
ind <- summary(fem)
time <- summary(fem_time)
tways <- summary(fem_twoways)
drsq <- data.frame(Individu=ind$r.squared,Time=time$r.squared,Twoways=tways$r.squared)

# Perbandingan
compare <- t(rbind(dsse,daic,dbic,dmape,drsq))
colnames(compare) <- c("SSE","AIC","BIC", "MAPE","R-Squared","Adj R-Squared")
compare

# FEM VS CEM
pooltest(cem, fem_twoways) # Keputusan pilih FEM

# REM dengan Generalized Least Square
rem_gls <- plm(PK~IPM+RPP+TPT+UM+GR, data = Dataset, 
               index = c("Kota", "Tahun"), 
               effect = "twoways", model = "random", random.method = "nerlove")
summary(rem_gls)

#efek individu
plmtest(rem_gls,type = "bp", effect="individu")

#efek waktu 
plmtest(rem_gls,type = "bp", effect="time")

#efek twoways 
plmtest(rem_gls,type = "bp", effect="twoways")

tidy_ranef_ind <- tidy(ranef(rem_gls, effect="individual"))
colnames(tidy_ranef_ind) <- c("Provinsi", "Pengaruh Acak Individu")
tidy_ranef_ind

tidy_ranef_time <- tidy(ranef(rem_gls, effect="time"))
colnames(tidy_ranef_time) <- c("Tahun", "Pengaruh Acak Waktu")
tidy_ranef_time

# FEM VS REM
# Uji Haussman
phtest(fem_twoways, rem_gls) # Pilih model FEM

# Uji diagnostik residu
# Uji normalitas
ks.test(fem_twoways$residuals, "pnorm", 
        mean=mean(fem_twoways$residuals), 
        sd=sd(fem_twoways$residuals))

shapiro.test(fem_twoways$residuals)

# Histogram
ggplot(as.data.frame(fem_twoways$residuals), aes(x = fem_twoways$residuals)) +
  geom_histogram(aes(y = after_stat(density)), color = "white", fill = "steelblue") +
  geom_density(color = "red", linewidth = 1) +
  theme_minimal()

# Uji Autokorelasi
pbgtest(fem_twoways)

# Uji Heteroskedastisitas
bptest(fem_twoways)

# Modeling GWPR
# Menentukan fungsi pembobot spasial terbaik
#adaptive bisquare
bw <- bw.GWPR(formula = PK~IPM+RPP+TPT+UM+GR, data = Dataset,
              index = c("ID", "Tahun"), SDF = Dataset.sdf, adaptive = TRUE,
              effect = "twoways", model = "within",
              kernel = "bisquare", longlat = FALSE,approach = "CV")

result <- GWPR(bw = bw, formula = PK~IPM+RPP+TPT+UM+GR, data = Dataset,
               index = c("ID", "Tahun"), SDF = Dataset.sdf, adaptive = TRUE,
               effect = "twoways", model = "within",
               kernel = "bisquare", longlat = FALSE)

# Define the Exponential kernel function
x = coordinates(as(Jatim_map,"Spatial"))[,1]
y = coordinates(as(Jatim_map,"Spatial"))[,2]
coords <-cbind(x,y)
jarak <-as.matrix(1/dist(coords))

bisquare_kernel <- function(dists, bw) {
  ifelse(dists < bw, (1-(dists/bw)^2)^2, 0)
}
weights <- bisquare_kernel(jarak, bw)
writexl::write_xlsx(data.frame(weights),"Pembobot Kernel Bisquare.xlsx")

# Perbandingan hasil
Compare2 <- data.frame(R.Squared = c(result$R2,summary(fem_twoways)$r.squared))
rownames(Compare2) <- c("GWPR","FEM","Adjusted R FEM")
Compare2

#Pembuatan jarak euclidean
n <- max(Dataset$ID) #jumlah wilayah
U <- Dataset$Long #data Longitude
V <- Dataset$Lat #data Latitude
d <- matrix(0,n,n)
for (i in 1:n) {
  for (j in 1:n) {
    d[i,j] <- sqrt(((U[i]-U[j])^2)+((U[i]-U[j])^2))
  }
}
d
writexl::write_xlsx(data.frame(d),"Hasil_Jarak_Euclidean.xlsx")

#----------------------------------------------------------------------------
#Export Hasil analisis GWPR ke Excel
GWPR.Result = st_as_sf(result$SDF)
writexl::write_xlsx(GWPR.Result,"hasil GWR Panel Terbaik.xlsx")

#---------------------------------------------------------------#
#    koefisien variabel (variabel IPM)
#---------------------------------------------------------------#
ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =IPM)) +
  scale_fill_gradient2(
    midpoint = 0, low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  ggtitle("Koefisien IPM")

#---------------------------------------------------------------#
#    koefisien variabel (variabel RPP)
#---------------------------------------------------------------#
ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =RPP)) +
  scale_fill_gradient2(
    midpoint = 0, low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  ggtitle("Koefisien RPP")

#---------------------------------------------------------------#
#    koefisien variabel (variabel TPT)
#---------------------------------------------------------------#
ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =TPT)) +
  scale_fill_gradient2(
    midpoint = 0, low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  ggtitle("Koefisien TPT")

#---------------------------------------------------------------#
#    koefisien variabel (variabel UM)
#---------------------------------------------------------------#
ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =UM)) +
  scale_fill_gradient2(
    midpoint = 0, low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  ggtitle("Koefisien UM")

#---------------------------------------------------------------#
#    koefisien variabel (variabel GR)
#---------------------------------------------------------------#
ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =GR)) +
  scale_fill_gradient2(
    midpoint = 0, low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  ggtitle("Koefisien GR")

#---------------------------------------------------------------#
#    Local R2
#---------------------------------------------------------------#
ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =Local_R2)) +
  scale_fill_gradient2(
    midpoint = mean(GWPR.Result$Local_R2), low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  ggtitle("Local R2")

#---------------------------------------------------------------#
#    signfikansi variabel (variabel IPM)
#---------------------------------------------------------------#
GWPR.Result$signfikansi_IPM<- NA
# Signifikan
GWPR.Result[(GWPR.Result$IPM_TVa <= -1.96 | GWPR.Result$IPM_TVa >= 1.96), "signfikansi_IPM"] <- "Signifikan"

# Tidak Signifikan
GWPR.Result[(GWPR.Result$IPM_TVa > -1.96 & GWPR.Result$IPM_TVa < 1.96), "signfikansi_IPM"] <- "Tidak Signifikan"

ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =signfikansi_IPM)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+ggtitle("Signifikansi IPM")

#---------------------------------------------------------------#
#    signfikansi variabel (variabel RPP)
#---------------------------------------------------------------#
GWPR.Result$signfikansi_RPP<- NA
# Signifikan
GWPR.Result[(GWPR.Result$RPP_TVa <= -1.96 | GWPR.Result$RPP_TVa >= 1.96), "signfikansi_RPP"] <- "Signifikan"

# Tidak Signifikan
GWPR.Result[(GWPR.Result$RPP_TVa > -1.96 & GWPR.Result$RPP_TVa < 1.96), "signfikansi_RPP"] <- "Tidak Signifikan"

ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =signfikansi_RPP)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+ggtitle("Signifikansi RPP")

#---------------------------------------------------------------#
#    signfikansi variabel (misal variabel TPT)
#---------------------------------------------------------------#
GWPR.Result$signfikansi_TPT <- NA
# Signifikan
GWPR.Result[(GWPR.Result$TPT_TVa <= -1.96 | GWPR.Result$TPT_TVa >= 1.96), "signfikansi_TPT"] <- "Signifikan"

# Tidak Signifikan
GWPR.Result[(GWPR.Result$TPT_TVa > -1.96 & GWPR.Result$TPT_TVa < 1.96), "signfikansi_TPT"] <- "Tidak Signifikan"

ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =signfikansi_TPT)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+ggtitle("Signifikansi TPT")

#---------------------------------------------------------------#
#    signfikansi variabel (misal variabel UM)
#---------------------------------------------------------------#
GWPR.Result$signfikansi_UM <- NA
# Signifikan
GWPR.Result[(GWPR.Result$UM_TVa <= -1.96 | GWPR.Result$UM_TVa >= 1.96), "signfikansi_UM"] <- "Signifikan"

# Tidak Signifikan
GWPR.Result[(GWPR.Result$UM_TVa > -1.96 & GWPR.Result$UM_TVa < 1.96), "signfikansi_UM"] <- "Tidak Signifikan"

ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =signfikansi_UM)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+ggtitle("Signifikansi UM")

#---------------------------------------------------------------#
#    signfikansi variabel (misal variabel GR)
#---------------------------------------------------------------#
GWPR.Result$signfikansi_GR <- NA
# Signifikan
GWPR.Result[(GWPR.Result$GR_TVa <= -1.96 | GWPR.Result$GR_TVa >= 1.96), "signfikansi_GR"] <- "Signifikan"

# Tidak Signifikan
GWPR.Result[(GWPR.Result$GR_TVa > -1.96 & GWPR.Result$GR_TVa < 1.96), "signfikansi_GR"] <- "Tidak Signifikan"

ggplot(data=GWPR.Result) +
  geom_sf(mapping=aes(fill =signfikansi_GR)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+ggtitle("Signifikansi GR")
