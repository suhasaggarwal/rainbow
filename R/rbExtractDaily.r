#' Extract daily values of weather variables from a PRISM raster data set
#' 
#' These functions extract values from interpolated weather rasters from the Parameter Regression on Independent Slopes (PRISM) data product (Daly et al. 2002). Depending on the set, PRISM rasters represent interpolated values of daily, monthly, or annual weather. Extractions are done at points (versus polygons or lines, for example). \cr\cr
#' The basic input is data frame or data frame-like object, with each row representing a location and having a particular date. The function also needs to be pointed toward a folder with PRISM data. The function assumes that the rasters are arranged in the manner in which they are normally supplied by the PRISM Climate Group. That is, for example, each daily/monthly/annual raster is inside a folder corresponding to the year, which in turn is in a folder named either \code{daily} or \code{monthly}, which is inside a folder named after the variable the rasters represent (e.g., \code{tmin}, \code{vpdmax}, etc.), which is inside a folder named either \code{an81} (daily values) or \code{lt81} (monthly and annual values). Not all years or variables are necessarily available, as a license must be purchased and data acquired from the PRISM Climate Group. \cr\cr
#' The function can extract values corresponding to the day/month/year of each record, plus (optionally) a user-specified window of time prior to the day/month/year of each record. For example, you could use this to extract daily climate data for a site collected on April 22, 2014, and all days prior up to 10 years (April 23, 2004). This function is really a fancy wrapper for \code{\link[terra]{extract}}, but it does the tough job of getting the directory structures right, pulling all needed rasters, and efficiently grouping records to speed extraction. \cr\cr
#' The function does not assume that data for all PRISM years are available (i.e., if you only have rights to the data from, say, 2015 to 2018), but it does assume that all relevant rasters for a particular year are available within each yearly folder. If rasters preceding a date only partially cover the window, then values for the part covered will be extracted. For example if you try to extract annual values for a window spanning 2010 to 2020 but only have available rasters for 1981 to 2018, then values for 2010 to 2018 will be extracted.
#'
#' @param prismDir Character. Path of the base directory in which the rasters are stored. Folder contains a subfolder named \code{an81} or \code{lt81}.
#' @param prismSet Character. Name of PRISM data subset. This can be either \code{an81} (daily weather data) or \code{lt81} (monthly and annual weather data).
#' @param rastVars Name of variable(s) to extract. Valid values are:
#' \itemize{
#'		\item	\code{ppt}: Precipitation.
#'		\item	\code{tmax}: Maximum temperature.
#'		\item	\code{tmin}: Minimum temperature.
#'		\item	\code{tmean}: Mean temperature.
#'		\item	\code{tdmean}: Dew-point temperature.
#'		\item	\code{vpdmax}: Maximum vapor-pressure.
#'		\item	\code{vpdmin}: Minimum vapor-pressure.
#' }
#' @param x Data frame.
#' @param yearField,monthField,domField Character. Names of columns in \code{x} with year, month, and day of month of observation. These must all be integers or \code{NA} values. You can specify these values or use argument \code{dateField} to specify a column with dates.
#' @param dateField Name of column in \code{x} with date of each record. Values must be of in YYYY-MM-DD (year-month-day of month) format \emph{or} already of class \code{\link{Date}}. See \code{\link[lubridate]{ymd}} or related functions for more help.
#' @param longLat Character vector with two elements. Names of columns in \code{x} with longitude and latitude (in that order). Coordinates will be assumed to be in the same coordinate reference system as the PRISM rasters.
#' @param windowYears,windowMonths,windowDays Integers >= 0. Extract data for this many years, months, and/or days before the day of each observation, \emph{plus} the day of observation. Note:
#' \itemize{
#'		\item If \code{prismSet} is \code{an81} (daily), then only \code{windowYears} and \code{windowDays} are used. Note that the difference between using, say, \code{windowYears = 1} and \code{windowDays = 365} is that the former can accommodate leap days whereas the latter just extracts the 365 days prior (which may be less than a full calendar year if the timespan encompasses a leap day).
#'		\item If \code{prismSet} is \code{lt81} (monthly), then only \code{windowYears} and \code{windowMonths} are used.
#'		\item If \code{prismSet} is \code{lt81} (annual), then only \code{windowYears}is used.
#' }
#' To get only data for the day/month/year of each record, set all of the respective \code{window} arguments to 0.
#' @param verbose Logical. If \code{TRUE} (default), show progress.
#'
#' @return Matrix with one row per row in \code{x}. \code{NA} values represent days/months/years that did not fall within the specified window or for which rasters were unavailable. 
#' @references Daly, C., Gibson, W.P., Taylor, G.H., Johnson, G.L., and Pasteris, P.  2002.  A knowledge-based approach to the statistical mapping of climate. \emph{Climate Research} 22:99-113. \href{http://dx.doi.org/10.3354/cr022099}{DOI: 10.3354/cr022099}
#' @export
#' @examples
#' \dontrun{
#' x <- data.frame(
#' 	long=rep(-97.66, 5),
#' 	lat=rep(38.37, 5),
#' 	date=c('2015-12-31', '1981-01-05', '2020-12-01',
#' 	'2019-01-05', '1895-05-01')
#' )
#' 
#' y <- rbExtractDaily(
#' 	prismDir = 'F:/ecology/Climate/PRISM acquired in 2020',
#' 	prismSet = 'an81',
#' 	rastVars = 'tmin',
#' 	x = x,
#' 	yearField = NULL,
#' 	monthField = NULL,
#' 	domField = NULL,
#' 	dateField = 'date',
#' 	longLat = c('long', 'lat'),
#' 	windowYears = 0,
#' 	windowDays = 7,
#' 	verbose = TRUE
#' )
#' 
#' y <- rbExtractMonthly(
#' 	prismDir = 'F:/ecology/Climate/PRISM acquired in 2020',
#' 	prismSet = 'lt81',
#' 	rastVars = 'tmin',
#' 	x = x,
#' 	yearField = NULL,
#' 	monthField = NULL,
#' 	domField = NULL,
#' 	dateField = 'date',
#' 	longLat = c('long', 'lat'),
#' 	windowYears = 1,
#' 	windowMonths = 1,
#' 	verbose = TRUE
#' )
#' 
#' y <- rbExtractYearly(
#' 	prismDir = 'F:/ecology/Climate/PRISM acquired in 2020',
#' 	prismSet = 'lt81',
#' 	rastVars = 'tmin',
#' 	x = x,
#' 	yearField = NULL,
#' 	monthField = NULL,
#' 	domField = NULL,
#' 	dateField = 'date',
#' 	longLat = c('long', 'lat'),
#' 	windowYears = 3,
#' 	verbose = TRUE
#' )
#' 
#' }
#' @export

rbExtractDaily <- function(
	prismDir,
	prismSet = 'an81',
	rastVars,
	x,
	yearField = 'year',
	monthField = 'month',
	domField = 'day',
	dateField = NULL,
	longLat = NULL,
	windowYears = 0,
	windowDays = 7,
	verbose = TRUE
) {

	### get dates
	dates <- .getDates(x, yearField, monthField, domField, dateField)
	
	### DAILY PRISM extraction
	##########################
	if (tolower(prismSet) == 'an81') {

		recordDates <- dates
		dates <- unique(sort(dates))

		### by VARIABLE
		for (rastVar in rastVars) {
		
			### get available years
			yearsAvail <- list.dirs(paste0(prismDir, '/', prismSet, '/', rastVar, '/daily'), full.names=FALSE, recursive=FALSE)
			yearsAvail <- as.integer(yearsAvail)

			earliestRasterDate <- lubridate::make_date(min(yearsAvail), '01', '01')
			latestRasterDate <- lubridate::make_date(max(yearsAvail), '12', '31')

			### create data frame to store extracted values
			daysNeeded <- windowDays + ceiling(365.25 * windowYears)
			thisOut <- matrix(NA, nrow=nrow(x), ncol=daysNeeded + 1)
			colnames(thisOut) <- paste0(prismSet, '_', rastVar, '_', daysNeeded:0, 'daysPrior')

			### extract by date
			for (countDate in seq_along(dates)) {
			
				date <- dates[countDate]
				if (verbose) print(date); flush.console()
			
				# get date of window start
				windowStartDate <- date
				lubridate::year(windowStartDate) <- lubridate::year(windowStartDate) - windowYears
				if (windowYears > 0) lubridate::day(windowStartDate) <- lubridate::day(windowStartDate) + 1
				lubridate::day(windowStartDate) <- lubridate::day(windowStartDate) - windowDays
				windowStartDate <- max(windowStartDate, earliestRasterDate)
			
				windowEndDate <- date
				windowEndDate <- min(windowEndDate, latestRasterDate)
			
				# date falls within PRISM data set range
				if (date >= earliestRasterDate & windowStartDate <= latestRasterDate) {
			
					datesNeeded <- seq(as.POSIXct(windowStartDate), as.POSIXct(windowEndDate), by='1 days')
					datesNeeded <- as.Date(datesNeeded)
			
					yearsNeeded <- lubridate::year(datesNeeded)
					monthsNeeded <- lubridate::month(datesNeeded)
					daysNeeded <- lubridate::day(datesNeeded)
			
					rastDates <- paste0(yearsNeeded, sprintf('%02.0f', monthsNeeded), sprintf('%02.0f', daysNeeded)) 
			
					rastsNeeded <- paste0(prismDir, '/', prismSet, '/', rastVar, '/daily/', yearsNeeded, '/prism_', rastVar, '_us_30s_', rastDates, '.tif')

					# get rasters
					rasts <- terra::rast(rastsNeeded)
				
					# extract to records with this date
					index <- which(recordDates == date)
					locs <- x[index, longLat]
					ext <- terra::extract(rasts, locs)

					ext <- ext[ , 2:ncol(ext), drop=FALSE]

					# remember
					lastCol <- ncol(thisOut) - as.integer(date - windowEndDate)
					firstCol <- lastCol - ncol(ext) + 1
					thisOut[index, firstCol:lastCol] <- ext
			
				} # if record date falls within years available in PRISM
			
			} # next date

			out <- if (exists('out', inherits=FALSE)) {
				cbind(out, thisOut)
			} else {
				thisOut
			}
		
		} # next variable
		
	}
	
	while (all(is.na(out[ , 1]))) out <- out[ , 2:ncol(out), drop=FALSE]
	out

}

#' @describeIn rbExtractDaily Extract monthly values of weather variables from a PRISM raster data set
#' @export
rbExtractMonthly <- function(
	prismDir,
	prismSet = 'lt81',
	rastVars,
	x,
	yearField = 'year',
	monthField = 'month',
	domField = 'day',
	dateField = NULL,
	longLat = NULL,
	windowYears = 0,
	windowMonths = 0,
	verbose = TRUE
) {

	### get dates
	dates <- .getDates(x, yearField, monthField, domField, dateField)

	### allocate window monthly "years" to years
	windowYears <- windowYears + floor(windowMonths / 12)
	windowMonths <- windowMonths %% 12

	### DAILY PRISM extraction
	##########################
	if (tolower(prismSet) == 'lt81') {

		recordYears <- lubridate::year(dates)
		recordMonths <- lubridate::month(dates)
		recordYearsMonths <- cbind(recordYears, recordMonths)

		uniqueYearsMonths <- recordYearsMonths[!duplicated(recordYearsMonths), ]
		uniqueYearsMonths <- uniqueYearsMonths[order(uniqueYearsMonths[ , 2]), ]
		uniqueYearsMonths <- uniqueYearsMonths[order(uniqueYearsMonths[ , 1]), ]
		
		### by VARIABLE
		for (rastVar in rastVars) {
		
			### get available years
			yearsAvail <- list.dirs(paste0(prismDir, '/', prismSet, '/', rastVar, '/monthly'), full.names=FALSE, recursive=FALSE)
			yearsAvail <- as.integer(yearsAvail)

			earliestRastYearMonth <- cbind(min(yearsAvail), 1)
			latestRastYearMonth <- cbind(max(yearsAvail), 12)

			### create data frame to store extracted values
			monthsNeeded <- windowMonths + (12 * windowYears)
			thisOut <- matrix(NA, nrow=nrow(x), ncol=monthsNeeded + 1)
			colnames(thisOut) <- paste0(prismSet, '_', rastVar, '_', monthsNeeded:0, 'monthsPrior')

			### extract by date
			for (i in 1:nrow(uniqueYearsMonths)) {
			
				yearMonth <- uniqueYearsMonths[i, , drop=FALSE]
				if (verbose) print(paste0(yearMonth[1, 1], '-', sprintf('%02.0f', yearMonth[1, 2]))); flush.console()
			
				# get date of window start
				windowStart <- yearMonth
				windowStart <- if (windowMonths < windowStart[1, 2]) {
					windowStart[1, ] - cbind(windowYears, windowMonths)
				} else {
					cbind(windowStart[1, 1] - windowYears - 1, 12 - (windowMonths - windowStart[1, 2]))
				}
				
				if (windowStart[1, 1] < earliestRastYearMonth[1, 1]) windowStart <- earliestRastYearMonth
				
				windowEnd <- yearMonth
				if (windowEnd[1, 1] > latestRastYearMonth[1, 1]) windowEnd <- latestRastYearMonth
			
				# date falls within PRISM data set range
				if (windowStart[1, 1] >= earliestRastYearMonth[1, 1] & windowStart[1, 1] <= latestRastYearMonth[1, 1]) {
			
					# get years and months needed
					yearsMonthsNeeded <- windowStart
					n <- 1
					while (any(yearsMonthsNeeded[n, ] < windowEnd)) {
					
						deltaMonth <- yearsMonthsNeeded[n, 2] + 1
						newYear <- if (deltaMonth == 13) { yearsMonthsNeeded[n, 1] + 1 } else { yearsMonthsNeeded[n, 1] }
						newMonth <- if (deltaMonth == 13) { 1 } else { yearsMonthsNeeded[n, 2] + 1 }
						yearsMonthsNeeded <- rbind(yearsMonthsNeeded, cbind(newYear, newMonth))
						n <- n + 1
					
					}
					
					yearsNeeded <- yearsMonthsNeeded[ , 1]
					monthsNeeded <- yearsMonthsNeeded[ , 2]
					
					rastDates <- paste0(yearsNeeded, sprintf('%02.0f', monthsNeeded)) 
					rastsNeeded <- paste0(prismDir, '/', prismSet, '/', rastVar, '/monthly/', yearsNeeded, '/prism_', rastVar, '_us_30s_', rastDates, '.tif')

					# get rasters
					rasts <- terra::rast(rastsNeeded)
				
					# extract to records with this date
					index <- which(yearMonth[1, 1] == recordYearsMonths[ , 1] & yearMonth[1, 2] == recordYearsMonths[ , 2])
					locs <- x[index, longLat]
					ext <- terra::extract(rasts, locs)

					ext <- ext[ , 2:ncol(ext), drop=FALSE]

					# remember
					deltaMonths <- sum(yearMonth * cbind(12, 1)) - sum(windowEnd * cbind(12, 1))
					lastCol <- ncol(thisOut) - deltaMonths
					firstCol <- lastCol - ncol(ext) + 1
					thisOut[index, firstCol:lastCol] <- ext
			
				} # if record date falls within years available in PRISM
			
			} # next date

			out <- if (exists('out', inherits=FALSE)) {
				cbind(out, thisOut)
			} else {
				thisOut
			}
		
		} # next variable
		
	}
	
	while (all(is.na(out[ , 1]))) out <- out[ , 2:ncol(out), drop=FALSE]
	out

}

#' @describeIn rbExtractDaily Extract annual values of weather variables from a PRISM raster data set
#' @export
rbExtractYearly <- function(
	prismDir,
	prismSet = 'lt81',
	rastVars,
	x,
	yearField = 'year',
	monthField = 'month',
	domField = 'day',
	dateField = NULL,
	longLat = NULL,
	windowYears = 0,
	verbose = TRUE
) {

	### get dates
	dates <- .getDates(x, yearField, monthField, domField, dateField)

	### extract annual data
	if (tolower(prismSet) == 'lt81') {

		# get record years
		recordYears <- lubridate::year(dates)
		uniqueYears <- sort(unique(recordYears))

		### by VARIABLE
		for (rastVar in rastVars) {
		
			### get available years
			yearsAvail <- list.dirs(paste0(prismDir, '/', prismSet, '/', rastVar, '/monthly'), full.names=FALSE, recursive=FALSE)
			yearsAvail <- as.integer(yearsAvail)
			
			earliestRasterYear <- min(yearsAvail)
			latestRasterYear <- max(yearsAvail)

			### create matrix to store extracted values
			thisOut <- matrix(NA, nrow=nrow(x), ncol=windowYears + 1)
			colnames(thisOut) <- paste0(prismSet, '_', rastVar, '_', windowYears:0, 'yearsPrior')

			### extract by date
			for (i in seq_along(uniqueYears)) {
			
				year <- uniqueYears[i]
				if (verbose) cat(' ', as.integer(year), '\n'); flush.console()
			
				# get date of window start
				startYear <- year - windowYears
				
				# date falls within PRISM data set range
				if (year >= earliestRasterYear & startYear <= latestRasterYear) {
			
					# get years of extraction
					windowStartYear <- max(earliestRasterYear, startYear)
					windowEndYear <- min(latestRasterYear, year)
					
					yearsNeeded <- windowStartYear:windowEndYear
					
					rastsNeeded <- paste0(prismDir, '/', prismSet, '/', rastVar, '/monthly/', yearsNeeded, '/prism_', rastVar, '_us_30s_', yearsNeeded, '.tif')

					# get rasters
					rasts <- terra::rast(rastsNeeded)
				
					# extract to records with this date
					index <- which(recordYears == year)
					locs <- x[index, longLat]
					ext <- terra::extract(rasts, locs)

					ext <- ext[ , 2:ncol(ext), drop=FALSE]

					# remember
					deltaYears <- if (year < latestRasterYear) {
						0
					} else {
						year - latestRasterYear
					}
					
					lastCol <- ncol(thisOut) - deltaYears
					firstCol <- lastCol - ncol(ext) + 1
					thisOut[index, firstCol:lastCol] <- ext
			
				} # if record date falls within years available in PRISM
			
			} # next date

			out <- if (exists('out', inherits=FALSE)) {
				cbind(out, thisOut)
			} else {
				thisOut
			}
		
		} # next variable

	} # if MONTHLY PRISM

	while (all(is.na(out[ , 1]))) out <- out[ , 2:ncol(out), drop=FALSE]
	out

}

### get dates from input data frame
.getDates <- function(x, yearField, monthField, domField, dateField) {
	
	if (!is.null(yearField) & !is.null(monthField) & !is.null(domField)) {
		dates <- lubridate::make_date(
			year=x[ , yearField, drop=TRUE],
			month=x[ , monthField, drop=TRUE],
			day=x[ , domField]
		)
	} else if (!is.null(dateField)) {
		if (class(x[ , dateField]) != 'Date') {
			dates <- lubridate::ymd(x[ , dateField, drop=TRUE])
		} else {
			dates <- x[ , dateField, drop=TRUE]
		}
	} else {
		stop('You must specify either "yearField", "monthField", and "domField" OR "dateField."')
	}
	
	dates
	
}
