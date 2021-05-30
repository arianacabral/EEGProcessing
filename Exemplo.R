# Indicar o diretório de trabalho 

setwd("D:/Ariana/Documentos/2021/PSB/EXEMPLO/EEGProcessing/")

# -------------------------------------------------------------------------

source("EEGProcessing.r")

filename <- file.choose() # escolha o arquivo desejado

EEG.Processing(filename) # processamento dos sinais EEG