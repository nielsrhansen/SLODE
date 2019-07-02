##-----------------------------------------##
## Defining the system of Hynne et al 2001 ##
##-----------------------------------------##

## Initial states 
x0 <- c(
  Glc = 0.573074, 
  Glcx = 1.55307, 
  G6P = 4.2,
  F6P = 0.49,
  FBP = 4.64,
  GAP = 0.115,
  DHAP = 2.95,
  BPG = 0.00027,
  PEP = 0.04,
  Pyr = 8.7,
  ACA = 1.48153,
  EtOH = 19.2379,
  EtOHx = 16.4514,
  Glyc = 4.196,
  Glycx = 1.68478,
  ACAx = 1.28836,
  CNx = 5.20358,
  ATP = 2.1,
  ADP = 1.5,
  AMP = 0.33,
  NADH = 0.33,
  NADp = 0.65)


## Stoichiometric matrix
CC <- matrix(0, nrow = 24, ncol = length(x0))
colnames(CC) <- names(x0)
rownames(CC) <- c("inGlc", "GlcTrans", "HK", "PGI", "PFK", "ALD", "TIM", "GAPDH", "lpPEP", "PK", "PDC", "ADH", "difEtOH", "outEtOH", "lpGlyc", "difGlyc", "outGlyc", "difACA", "outACA", "lacto", "inCN", "storage", "consum", "AK")

CC["inGlc", "Glcx"] <- 1
CC["GlcTrans", c("Glcx")] <- -1;  CC["GlcTrans", c("Glc")] <- 1
CC["HK", c("Glc", "ATP")] <- -1;  CC["HK", c("G6P", "ADP")] <- 1
CC["PGI", c("G6P")] <- -1;        CC["PGI", c("F6P")] <- 1
CC["PFK", c("G6P")] <- -1;        CC["PFK", c("F6P")] <- 1
CC["ALD", c("FBP")] <- -1;        CC["ALD", c("GAP", "DHAP")] <- 1
CC["TIM", c("DHAP")] <- -1;       CC["TIM", c("GAP")] <- 1
CC["GAPDH", c("GAP", "NADp")] <- -1;  CC["GAPDH", c("BPG", "NADH")] <- 1
CC["lpPEP", c("BPG", "ADP")] <- -1;   CC["lpPEP", c("PEP", "ATP")] <- 1
CC["PK", c("PEP", "ADP")] <- -1;      CC["PK", c("Pyr", "ATP")] <- 1
CC["PDC", c("Pyr")] <- -1;            CC["PDC", c("ACA")] <- 1
CC["ADH", c("ACA", "NADH")] <- -1;    CC["ADH", c("EtOH", "NADp")] <- 1
CC["difEtOH", c("EtOH")] <- -1;       CC["difEtOH", c("EtOHx")] <- 1
CC["outEtOH", c("EtOHx")] <- -1;
CC["lpGlyc", c("DHAP", "NADH")] <- -1;  CC["lpGlyc", c("Glyc", "NADp")] <- 1
CC["difGlyc", c("Glyc")] <- -1;         CC["difGlyc", c("Glycx")] <- 1
CC["outGlyc", c("Glycx")] <- -1;
CC["difACA", c("ACA")] <- -1;           CC["difACA", c("ACAx")] <- 1
CC["outACA", c("ACAx")] <- -1;
CC["lacto", c("ACAx", "CNx")] <- -1;
CC["inCN", c("CNx")] <- -1;
CC["storage", c("G6P", "ATP")] <- -1;   CC["storage", c("ADP")] <- 1
CC["consum", c("ATP")] <- -1;           CC["consum", c("ADP")] <- 1
CC["AK", c("ATP", "AMP")] <- -1;        CC["AK", c("ADP")] <- 2



## Constants (forward/reverse)
kf <- c("0" = 0.048, "9" = 4.43866e5, "13" = 1.67200e1, "16" = 1.9, "18" = 2.47e1, "20" = 2.83828e-3, 
  "22" = 2.25932e0, "23" = 3.20760e0, "24" = 4.32900e2)
kr <- c("0" = 0.048, "9" = 1.52662e5, "13" = 1.67200e1, "16" = 1.9, "18" = 2.47e1, "20" = 2.83828e-3, 
  "22" = 2.25932e0, "23" = 3.20760e0, "24" = 1.33333e2)

Vf <- c("2" = 1.01496e3, "3" = 5.17547e1, "4" = 4.96042e2, "5" = 4.54327e1, "6" = 2.20782e3, "6f" = 5, "6r" = 5,
  "7" = 1.16365e2, "8" = 8.33858e2, "10" = 3.43096e2, "11" = 5.31328e1, "12" = 8.98023e1, "15" = 8.14797e1)
Vr <- c("2" = 1.01496e3, "3" = 5.17547e1, "4" = 4.96042e2, "5" = 4.54327e1, "6" = 1.10391e4, "6f" = 5, "6r" = 5,
  "7" = 1.16365e2, "8" = 8.33858e2, "10" = 3.43096e2, "11" = 5.31328e1, "12" = 8.98023e1, "15" = 8.14797e1)

K <- c("2Glc" = 1.7, "2IG6P" = 1.2, "2IIG6P" = 7.2, "P2" = 1,
  "3ATP" = 0.1, "3Glc" = 0, "3dGlc" = 0.37, 
  "4G6P" = 0.8, "4F6P" = 0.15, "4eq" = 0.021,
  "5I" = 0.021, "5II" = 0.15,
  "6eq" = 0.081, "6FBP" = 0.3, "6GAP" = 4.0, "6DHAP" = 2.0, "6IGAP" = 10.0,
  "7DHAP" = 1.23, "7GAP" = 1.27, "7eq" = 0.055,
  "8GAP" = 0.6, "8BPG" = 0.01, "8NAD" = 0.1, "8NADH" = 0.06, "8eq" = 0.06,
  "10ADP" = 0.17, "10PEP" = 0.2, 
  "11" = 0.3,
  "12ACA" = 0.71, "12NADH" = 0.1,
  "13" = 1.9 * 8.8,
  "15NADH" = 0.13,
  "15DHAP" = 25,
  "15INADH" = 0.034,
  "15INAD" = 0.13,
  "16" = 1.9,
  "18" = 1.9 * 13)

e <- c("yvol" = 59, "Glcx0" = 18.5, "CNx0" = 5.60)





## Complex names
complex_names <- c(
  "0", "Glcx", 
  "Glc", "G6P", "Glc,Glcx", "G6P,Glcx", "G6P,Glc", "G6P,Glc,Glcx",
  "ATP", "ATP,Glc",
  "F6P",
  "AMP2,F6P2", "AMP2", "ATP2",
  "FBP", "GAP", "DHAP", "FBP,GAP", "DHAP,GAP",
  
  "GAP,NADp", "BPG,NADH", "BPG", "NADp", "NADH", "GAP,NADH", "BPG,NADp",
  "ADP,BPG", "ATP,PEP",
  "ADP", "PEP", "ADP,PEP",
  "Pyr",
  "ACA,NADH", "ACA",
  "EtOH", "EtOHx",
  
  "DHAP,NADH", "DHAP,NADp",
  "Glyc", "Glycx",
  
  "ACAx",
  
  "ACAx,CNx",
  "CNx",
  "ATP,G6P",
  
  "AMP,ATP", "ADP2")

stopifnot(length(complex_names) == length(unique(complex_names)))

## Power matrix A
A <- matrix(0, nrow = length(complex_names), ncol = length(x0))
rownames(A) <- complex_names
colnames(A) <- names(x0)
for (i in 2:nrow(A)) {
  cn <- rownames(A)[i]                    # get complex name
  nn <- strsplit(cn, ",")[[1]]            # split by ,
  no <- substr(nn, nchar(nn), nchar(nn))  # get last char (may be number)
  has_no <- no %in% 1:10                  # is the last a number
  nn[has_no] <- substr(nn[has_no], 1, nchar(nn[has_no]) - 1)  # remove if so
  ex <- rep(1, length(no))                # exponent
  ex[has_no] <- as.numeric(no[has_no])    # exponent if it has number
  A[i, nn] <- ex
}



## Coefficient matrices, note if reaction goes both ways, the share denominator, so let K1 have negatives
K1 <- matrix(0, nrow = nrow(CC), ncol = nrow(A))
rownames(K1) <- rownames(CC); colnames(K1) <- rownames(A)
K2 <- K1

# order first -> (positive), then <- (negative) (if exists)
K1[1, c("0")] <- e["yvol"] * kf["0"] * e["Glcx0"];    
K1[1, c("Glcx")] <- -e["yvol"] * kf["0"];

K1[2, c("Glcx", "Glc,Glcx")] <- 
  c(Vf["2"] / K["2Glc"], Vf["2"] * K["P2"] / K["2Glc"]^2);    
K1[2, c("Glc", "Glc,Glcx")] <- K1[2, c("Glc", "Glc,Glcx")] - 
  c(Vf["2"] / K["2Glc"], Vr["2"] * K["P2"] / K["2Glc"]^2);
K2[2, c("0", "Glc", "Glcx", "G6P", "Glc,Glcx", "G6P,Glcx", "G6P,Glc", "G6P,Glc,Glcx")] <- 
  c(1, (1 + K["P2"]) / K["2Glc"], (1 + K["P2"]) / K["2Glc"], 1 / K["2IG6P"], 2 * K["P2"] / K["2Glc"]^2, 
    K["P2"] / (K["2Glc"] * K["2IG6P"]), 1 / (K["2Glc"] * K["2IIG6P"]), K["P2"] / (K["2Glc"]^2 * K["2IG6P"]))

K1[3, c("ATP,Glc")] <- c(Vf["3"])
K2[3, c("ATP", "Glc", "ATP,Glc")] <- c(K["3Glc"], K["3ATP"], 1)
K1[3, ] <- K1[3, ] / (K["3dGlc"] * K["3ATP"])
K2[3, ] <- K2[3, ] / (K["3dGlc"] * K["3ATP"])

K1[4, c("G6P")] <- c(Vf["4"])
K1[4, c("F6P")] <- K1[4, c("F6P")] - c(Vr["4"] / K["4eq"])
K2[4, c("G6P", "F6P")] <- c(1, K["4G6P"] / K["4F6P"])
K1[4, ] <- K1[4, ] / K["4G6P"]
K2[4, ] <- K2[4, ] / K["4G6P"]

eps <- 1e-5
K1[5, c("AMP2,F6P2")] <- c(Vf["5"])
K2[5, c("AMP2", "ATP2", "AMP2,F6P2")] <- c(K["5I"], K["5II"], 1)
K1[5, ] <- K1[5, ] / eps
K2[5, ] <- K2[5, ] / eps

K1[6, c("FBP")] <- c(Vf["6"])
K1[6, c("DHAP,GAP")] <- K1[6, c("DHAP,GAP")] - c(Vr["6"] / K["6eq"])
K2[6, c("FBP", "GAP", "DHAP", "FBP,GAP", "DHAP,GAP")] <- 
  c(1, K["6DHAP"] * Vf["6"] / (K["6eq"] * Vr["6"]), K["6GAP"] * Vf["6"] / (K["6eq"] * Vr["6"]),
    1 / K["6IGAP"], Vf["6"] / (K["6eq"] * Vr["6"]))
K1[6, ] <- K1[6, ] / K["6FBP"]
K2[6, ] <- K2[6, ] / K["6FBP"]

K1[7, c("DHAP")] <- c(Vf["7"])
K1[7, c("GAP")] <- K1[7, c("GAP")] - c(Vr["7"] / K["7eq"])
K2[7, c("DHAP", "GAP")] <- c(1, K["7DHAP"] * K["7GAP"])
K1[7, ] <- K1[7, ] / K["7DHAP"]
K2[7, ] <- K2[7, ] / K["7DHAP"]

K1[8, c("GAP,NADp")] <- c(Vf["8"] / (K["8GAP"] * K["8NAD"]))
K1[8, c("BPG,NADH")] <- K1[8, c("BPG,NADH")] - c(Vf["8"] / (K["8GAP"] * K["8NAD"] * K["8eq"]))
K2[8, c("GAP", "BPG", "NADp", "NADH", "GAP,NADp", "BPG,NADp", "GAP,NADH", "BPG,NADH")] <- 
  1 / c(K[c("8GAP", "8BPG", "8NAD", "8NADH")], K["8GAP"] * K["8NAD"], K["8BPG"] * K["8NAD"], K["8GAP"] * K["8NADH"], K["8BPG"] * K["8NADH"])

K1[9, c("ADP,BPG")] <- c(kf["9"])
K1[9, c("ATP,PEP")] <- K1[9, c("ATP,PEP")] - kr["9"]

K1[10, c("ADP,PEP")] <- c(Vf["10"])
K2[10, c("PEP", "ADP", "ADP,PEP")] <- c(K["10ADP"], K["10PEP"], 1)
K1[10, ] <- K1[10, ] / (K["10PEP"] * K["10ADP"])
K2[10, ] <- K2[10, ] / (K["10PEP"] * K["10ADP"])

K1[11, c("Pyr")] <- Vf["11"]
K2[11, c("Pyr")] <- 1
K1[11, ] <- K1[11, ] / K["11"]
K2[11, ] <- K2[11, ] / K["11"]

K1[12, c("ACA,NADH")] <- Vf["12"]
K2[12, c("NADH", "ACA", "ACA,NADH")] <- c(K["12ACA"], K["12NADH"], 1)
K1[12, ] <- K1[12, ] / (K["12ACA"] * K["12NADH"]) 
K2[12, ] <- K2[12, ] / (K["12ACA"] * K["12NADH"])

K1[13, c("EtOH")] <- c(kf["13"])
K1[13, c("EtOHx")] <- K1[13, c("EtOHx")] - kr["13"]

K1[14, c("EtOHx")] <- e["yvol"] * kf["0"]

K1[15, c("DHAP,NADH")] <- Vf["15"]
K2[15, c("NADH", "NADp", "DHAP,NADH", "DHAP", "DHAP,NADp")] <- 
  c(K["15DHAP"], K["15DHAP"] * K["15INADH"] / K["15INAD"], 1, K["15NADH"], K["15NADH"] / K["15INAD"])
K1[15, ] <- K1[15, ] / (K["15DHAP"] * K["15INADH"]) 
K2[15, ] <- K2[15, ] / (K["15DHAP"] * K["15INADH"]) 

K1[16, c("Glyc")] <- c(kf["16"])
K1[16, c("Glycx")] <- K1[16, c("Glycx")] - c(kr["16"])

K1[17, c("Glycx")] <- e["yvol"] * kf["0"]

K1[18, c("ACA")] <- c(kf["18"])
K1[18, c("ACAx")] <- K1[18, c("ACAx")] - c(kr["18"])

K1[19, c("ACAx")] <- e["yvol"] * kf["0"]

K1[20, c("ACAx,CNx")] <- e["yvol"] * kf["20"]

K1[21, c("CNx")] <- e["yvol"] * kf["0"]
K1[21, "0"] <- K1[21, "0"] -  e["yvol"] * kf["0"] * e["CNx0"]

K1[22, c("ATP,G6P")] <- kf["22"]

K1[23, c("ATP")] <- kf["23"]

K1[24, c("AMP,ATP")] <- c(kf["24"])
K1[24, c("ADP2")] <- K1[24, c("ADP2")] - c(kr["24"])



rm(cn, complex_names, e, eps, ex, has_no, i, K, kf, kr, nn, no, Vf, Vr)


# Standardise
K1 <- K1 / sqrt(rowSums(K1^2))
K2 <- K2 / sqrt(rowSums(K1^2))
