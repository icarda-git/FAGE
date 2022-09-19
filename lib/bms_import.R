#' Title line
#' 
#' Description paragraph
#' 
#' @param study.name  The name of the trial as listed in BMS (e.g., "IDYT42") 
#' @param trait.name  The trait name as defined in the BMS ontology (e.g., "GY_Calc_kgha")
#' @param bms.crop    The crop name in the BMS database (e.g., "wheat")
#' @param bms.program The name of the breeding program where trial dataset exists (e.g., "Wheat International Nurseries")
#' @param bms.user    The BMS username to login and get access to the backend database (optional, default is NULL). If value is NULL, then a login window will pop-up to insert the BMS username and password.
#' @param bms.pass    The BMS password to login and get access to the backend database (optional, default is NULL). If value is NULL, then a login window will pop-up to insert the BMS username and password.
#' @param bms.url     The URL of the BMS login page (default is "https://bms.icarda.org/ibpworkbench/")
#' 
#' @return description
#' 
#' @author 
#' Khaled Al-Shamaa. ICARDA
#' 
#' @examples 
#' # R code

bms_import <- function(study.name=NULL, trait.name=NULL, bms.crop=NULL, bms.program=NULL,
                       bms.url='https://bms.icarda.org/ibpworkbench/', bms.user=NULL, bms.pass=NULL){

  # config your BMS connection
  set_qbms_config(bms.url)
  
  # login using your BMS account (interactive mode)
  # you can pass BMS username and password as parameters (batch mode)
  login_bms(bms.user, bms.pass)
  
  # select wheat crop, you can list supported crops using list_crops() function
  set_crop(bms.crop)
  
  # select the breeding program by name 
  # you can list existing breeding programs in the selected crop using list_programs() function
  set_program(bms.program)
  
  # select a specific trial by name (e.g. IDYT37, IDYT38, IDYT39, IDYT40, IDYT41)
  # you can list all trials available in the selected breeding program using list_trials() function
  set_trial(study.name)
  
  # retrive multi-environment trial data
  trial.data <- get_trial_data()
  germplasm  <- get_germplasm_list()
  
  # get required data to upload Single-Site-Analysis BLUEs and summary stats files to BMS
  trial.data$TRIAL_INSTANCE <- gsub(".+ (\\d+)$", "\\1", trial.data$studyName)
  trial.data$studyId        <- debug_qbms()$state$trial_db_id
  #germplasm$GID

  # exclude the location if all trait values is NA
  #has.data   <- trial.data %>% group_by(studyName) %>% filter(!is.na((!!as.name(trait.name)))) %>% summarize(n())
  has.data   <- unique(na.omit(trial.data[,c('studyName',trait.name)])[,'studyName'])
  trial.data <- trial.data[trial.data$studyName %in% has.data,]
  
  # create a dummy environment name to handle the issue of multi trials in the same location
  trial.data <- trial.data %>% 
    arrange(locationName, studyName) %>% 
    group_by(locationName) %>% 
    mutate(envRank = rank(studyName, ties.method='min'), totalPlots = max(as.numeric(plotNumber))) %>%
    mutate(Env = ifelse(envRank==1,as.character(locationName), paste(locationName, 1+round(envRank/totalPlots))))
  
  # select required columns for further analysis steps 
  # Block data for incomplete block designs (e.g. Alpha design)
  # Row and Column data for spatial analysis
  trial.data <- trial.data[, c('studyId', 'TRIAL_INSTANCE', 'Env', 'germplasmName', 'entryType', 'replicate', 'blockNumber', 'X', 'Y', trait.name)]
  colnames(trial.data) <- c('studyId', 'TRIAL_INSTANCE', 'trial', 'gen', 'check', 'block', 'ibk', 'col', 'row', 'resp')
  
  # represent missing values of spatial structure properly
  trial.data$col[trial.data$col == 'null'] <- NA
  trial.data$row[trial.data$row == 'null'] <- NA
  
  # variable identifying checks (categories 0:TEST, 1:CHECK)
  trial.data$check <- as.character(trial.data$check)
  trial.data$check[trial.data$check == 'Test entry'] <- 0
  trial.data$check[trial.data$check == 'Check entry'] <- 1
  trial.data$check <- as.integer(trial.data$check)
  
  # tag checks in the germplasm list
  check.names <- unique(trial.data[trial.data$check == 1,]$gen)
  germplasm$check <- ifelse(germplasm$germplasmName %in% check.names, 1, 0)

  # drop national/local check (National Check)
  local.check <- c('nat check', 'national check', 'local check')
  trial.data  <- droplevels(subset(trial.data, !(tolower(trial.data$gen) %in% local.check)))
  germplasm   <- droplevels(subset(germplasm, !(tolower(germplasm$germplasmName) %in% local.check)))
  
  # save intermediate data for staging process
  data.file <- paste0(study.name, '_', trait.name, '_data.csv')
  geno.file <- paste0(study.name, '_germplasm.csv')
  
  # save dataset and germplasm list in intermediate files
  write.csv(as.data.frame(trial.data), file = data.file, row.names = FALSE)
  write.csv(germplasm, file = geno.file)
  
  return(list(data=data.file, germplasm=geno.file))
}
