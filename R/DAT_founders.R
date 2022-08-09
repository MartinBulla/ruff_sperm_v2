con = dbcon(db = "RUFFatSEEWIESEN")
x = dbq(con, "select * from FOUNDERS")
y = dbq(con, "select * from ADULTS where dead is not null")
xy = merge(x,y[,.(ID, date, dead)], by ='ID',all.x = TRUE)
xy[dead == 1,.(ID,arrival_date, date)]

z[, diff := as.numeric(difftime(date,arrival_date, unit = 'days'))]
z[diff<365]


summary(factor(x$population[x$arrival_date<as.POSIXct('2021-04-01')]))
48+49+55+55
summary(factor(x$population[x$arrival_date<as.POSIXct('2021-04-01') & is.na(x$comments)]))