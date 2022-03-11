#Jenny Smith

setwd(file.path(PROJHOME,"2017.11.07_MSLN_Correlation"))


tpm.1031 <- read.csv("~/RNA_seq_Analysis/0000.00.03_Expression_Matrices/TARGET_AML_AAML1031_dupGenesRemoved_TPM.csv", 
                     stringsAsFactors = FALSE, row.names = 1)

head(tpm.1031[,1:5])


orig <- read.csv("~/reference_mapping-files/AAML1031_TARGET_format_033118_v3_withPostInductionSurivalData.csv", 
                 stringsAsFactors = FALSE)


head(orig[,1:5])


# orig <- orig %>%
  # left_join(., dplyr::select(CDE.1031,USI,Reg.), by=c("Patient.registration.number"="Reg."))


addMSLN <- tpm.1031 %>% 
  rownames_to_column("Gene") %>%
  filter(Gene=="MSLN") %>%
  gather(USI,MSLN.in.TPM,-Gene) %>%
  mutate(Group=ifelse(grepl("^BM|^RO", USI), "NBM", "AML")) %>%
  left_join(., dplyr::select(CDE.1031,USI,Reg.), by=c("USI")) %>%
  dplyr::select(-Gene)

head(addMSLN)
dim(addMSLN)


addMSLN.NBM <- addMSLN %>%
  filter(is.na(Reg.))

dim(addMSLN.NBM)

addMSLN <- addMSLN %>%
  filter(! is.na(Reg.))

dim(addMSLN)


addMSLN.CDE <- orig %>%
  left_join(.,addMSLN, by=c("Patient.registration.number"="Reg.")) %>%
  bind_rows(., addMSLN.NBM) %>%
  dplyr::select(Patient.registration.number, USI,Group,MSLN.in.TPM, everything(),-Reg.)


  
dim(addMSLN.CDE)  
head(addMSLN.CDE)

write.csv(addMSLN.CDE, "~/AAML1031_TARGET_format_033118_v3_withPostInductionSurivalData_withMSLN.csv", row.names = FALSE)






