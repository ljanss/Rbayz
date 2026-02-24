# Download stock prices (standalone script)
#
# Download historical stock prices for one or more symbols from a data source
# (default: Yahoo Finance). Returns a named list of xts objects (one per symbol)
# unless a single symbol is supplied, in which case the xts object is returned.
#
# Usage (interactive):
#   prices <- download_stock_prices(c("AAPL", "MSFT"), from = "2020-01-01")
#
# Usage (script):
#   Rscript download_stock_prices.R AAPL,MSFT 2020-01-01 2020-12-31

if (!requireNamespace("quantmod", quietly = TRUE)) {
	stop("Package 'quantmod' is required. Install it with install.packages('quantmod').")
}

download_stock_prices <- function(symbols,
																	from = "2000-01-01",
																	to = Sys.Date(),
																	src = "yahoo",
																	auto.assign = FALSE) {
	if (missing(symbols) || length(symbols) == 0) {
		stop("symbols must be a non-empty character vector")
	}

	symbols <- as.character(symbols)

	result <- lapply(symbols, function(sym) {
		tryCatch(
			quantmod::getSymbols(
				sym,
				from = from,
				to = to,
				src = src,
				auto.assign = auto.assign
			),
			error = function(err) {
				stop("Failed to download data for '", sym, "': ", err$message)
			}
		)
	})

	names(result) <- symbols

	if (length(result) == 1) {
		return(result[[1]])
	}

	result
}

if (!interactive() && sys.nframe() == 0) {
	args <- commandArgs(trailingOnly = TRUE)
	if (length(args) < 1) {
		stop("Usage: Rscript download_stock_prices.R SYMBOLS[,SYMBOLS...] [FROM] [TO]")
	}
	symbols <- strsplit(args[1], ",")[[1]]
	from <- if (length(args) >= 2) args[2] else "2000-01-01"
	to <- if (length(args) >= 3) args[3] else Sys.Date()

	prices <- download_stock_prices(symbols, from = from, to = to)
	print(prices)
}
