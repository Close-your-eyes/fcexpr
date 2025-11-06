.book_ilab <- function(login_name = NULL,
                       login_pw = NULL,
                       link_week = "https://eu.ilabsolutions.com/schedules/270426#/schedule/week7/2023-09-22",
                       link_slot = "https://eu.ilabsolutions.com/schedules/270426#/schedule/week7/2023-09-22/events/new/2023-09-22%2018:00/2023-09-22%2019:00") {

  library(RSelenium)

  rD <- RSelenium::rsDriver(port = netstat::free_port(), browser = "firefox", chromever = NULL)
  remDr <- rD[["client"]]
  remDr$navigate("https://eu.ilabsolutions.com/account/login")

  Sys.sleep(5)
  username_box <- remDr$findElement("name", "login")
  username_box$sendKeysToElement(list(login_name))
  password_box <- remDr$findElement("name", "password")
  password_box$sendKeysToElement(list(login_pw, key = "enter"))

  Sys.sleep(8)
  remDr$navigate(link_week)

  curr_datetime <- lubridate::as_datetime(gsub(" CEST", "", Sys.time()))

  target_time <- gsub("\\%20", "", stringr::str_extract(link_slot, "\\%[:digit:]{1,}:00"))
  target_time_hour <- as.numeric(strsplit(target_time, ":")[[1]][1])
  target_date <- strsplit(strsplit(link_slot, "/new/")[[1]][2], "\\%")[[1]][1]

  target_datetime_to_book <- lubridate::as_datetime(paste0(Sys.Date(), " ", target_time, ":00"))
  message("Target time for placing the booking: ", target_datetime_to_book)
  sleeptime <- lubridate::seconds(lubridate::as.difftime(lubridate::interval(curr_datetime, target_datetime_to_book))) - 6
  message("Sleeptime until booking: ", lubridate::seconds_to_period(sleeptime))

  Sys.sleep(sleeptime) # ends 5 sec prior to target_time
  message(Sys.time(), " let's go.")

  # does it work for one digit hours?
  while(lubridate::hour(lubridate::as_datetime(Sys.time())) != target_time_hour) {
    # do nothing
  }

  # wait until sec 1? # no, not needed
  #while(round(lubridate::second(lubridate::as_datetime(Sys.time()))) != 1) {# donothing}


  remDr$navigate(link_slot)
  remDr$refresh()

  Sys.sleep(2) # worked
  buttons <- remDr$findElements("class", "button")
  buttons[[length(buttons)]]$clickElement()

  Sys.sleep(5) # worked
  save_button <- remDr$findElements("class", "positive")
  save_button[[length(save_button)]]$clickElement()



}
