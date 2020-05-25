# Required packages -------------------------------------------------------

# if (!require('pacman', quietly = T)) install.packages('pacman')
pacman::p_load('optparse', 'futile.logger', 'checkmate', 'pbapply', 'ggplot2', 'plotly', 'tidyr', 'plyr', 'dplyr', 'tibble', 'hues')


# CLI Parameters ----------------------------------------------------------

# List of parameters.
option_list <- list(
    optparse::make_option(c('-i', '--input'), action = 'store', default = NA, type = 'character', help = 'Folder containing the QC and logs of the RNA-Seq NextFlow run.'),
    optparse::make_option(c('-t', '--threads'), action = 'store', default = 1, type = 'integer', help = 'Number of threads.'),
    optparse::make_option(c('-v', '--verbose'), action = 'store_true', default = FALSE, help = 'Print additional progress messages.')
)

# Initialize optparse.
opt <- optparse::parse_args(
    optparse::OptionParser(
        option_list=option_list,
        description = 'Combines various QC and sequencing metrics from the RNA-Seq NextFlow pipeline.',
        epilogue = 'Please report any issues on https://github.com/UMCUGenetics/RNAseq-NF/issues/'
    )
)

# If verbose, print some additional progress messages.
if(opt$verbose){
    x <- futile.logger::flog.threshold(futile.logger::TRACE)
}else{
    x <- futile.logger::flog.threshold(futile.logger::INFO)
}

opt$input <- '~/test/testRNA/output/QC/'

# Helper functions --------------------------------------------------------

readSortMeRNA <- function(x){

    data <- readr::read_lines(x)

    # Retrieve the lines containing the number of reads passing (rRNA) and not passing (mRNA)
    data <- data[base::grepl('Total reads.*E-value threshold', data)] %>% base::trimws(.) %>% tibble::tibble(raw = .)

    # Clean-up data.
    data <- data %>%
        tidyr::separate(raw, c('type', 'value'), ' = ') %>%
        tidyr::separate(value, c('value', 'value.relative'), ' ') %>%
        dplyr::mutate(value.relative = base::gsub('[\\(\\)]', '', value.relative)) %>%
        dplyr::mutate(value = as.integer(value), value.relative = as.double(value.relative) / 100)

    data$sample <- base::gsub('_.*', '', base::basename(x))

    return(data)

}

readPreseq <- function(x){
    data <- readr::read_delim(x, delim = '\t', trim_ws = T, col_types = readr::cols())
    data$sample <- base::gsub('_.*', '', base::basename(x))

    return(data)
}

readTrimGalore <- function(x){

    data <- readr::read_lines(x)

    # Retrieve the lines summarizing the number of total reads and those containing adapters.
    data <- data[base::grepl('Total reads processed|Reads with adapters|Reads written|Total basepairs processed|Quality-trimmed|Total written|Number of sequence pairs removed', data)] %>% base::trimws(.) %>% tibble::tibble(raw = .)

    # Clean-up data.
    data <- data %>%
        tidyr::separate(raw, c('type', 'value'), ':') %>%
        tidyr::separate(value, c('value', 'value.relative'), '\\(', fill = 'right') %>%
        dplyr::mutate(value = base::gsub('[,bp]', '', value)) %>%
        dplyr::mutate_all(trimws) %>%
        dplyr::mutate(value.relative = base::gsub('[\\(\\)%]', '', value.relative)) %>%
        dplyr::mutate(value = as.integer(value), value.relative = as.double(value.relative) / 100) %>%
        dplyr::mutate(value.relative = ifelse(is.na(value.relative), 1, value.relative)) %>%
        dplyr::mutate(type = ifelse(grepl('shorter than the length', type), gsub('.*length cutoff', 'Below length cutoff', type), type))

    data$sample <- base::gsub('_.*', '', base::basename(x))
    data$file <- base::gsub('\\.fastq.*', '', base::basename(x))
    data$Direction <- ifelse(grepl('_R1_', data$file), 'R1', 'R2')

    return(data)

}

readTIN <- function(x){
    data <- readr::read_delim(x, delim = '\t', col_types = readr::cols()) %>%
        dplyr::mutate(sample = base::gsub('_.*', '', base::basename(Bam_file)))

    return(data)
}

readFlagstat.sambamba <- function(x){

    # Read the flagstat file.
    data <- tibble::tibble(value = readr::read_lines(x)) %>%
        dplyr::mutate(value = gsub('\\+.*', '', value)) %>%
        dplyr::mutate_all(trimws) %>% dplyr::mutate_all(as.integer) %>%
        dplyr::mutate(type = c('readsPassedQC', 'readsMultiMapped', 'readsSupplementary', 'readsDuplicated', 'readsMapped', 'readsPaired', 'readsFirstMate', 'readsSecondMate', 'readsProperPair', 'readsBothMatesMapped', 'readsSingleton', 'readsMateMappedOtherChr', 'readsMateMappedOtherChrWithQ5')) %>%
        dplyr::mutate(sample = base::gsub('_.*', '', base::basename(x)))

    # Add additional row showing unique reads.
    data <- base::rbind(data, base::data.frame(value = data[data$type == 'readsPassedQC',]$value - data[data$type == 'readsMultiMapped',]$value, type = 'Unique reads', sample = unique(data$sample)))

    # Return dataframe containing single flagstat file information.
    return(data)

}

# Check input parameters --------------------------------------------------

futile.logger::flog.debug('Checking input parameters')

checkmate::assertDirectory(opt$input, access = 'r')
checkmate::assertInteger(opt$threads)
checkmate::assertLogical(opt$verbose)


futile.logger::flog.info('Importing and parsing QC data.')


# TrimGalore --------------------------------------------------------------

futile.logger::flog.info('\t- TrimGalore')

TrimGalore <- list()
TrimGalore$files <- base::list.files(path = opt$input, full.names = T, recursive = T, pattern = '_trimming_report.txt$')
TrimGalore$data <- do.call(base::rbind, pbapply::pblapply(TrimGalore$files, readTrimGalore, cl = opt$threads))


# SortMeRNA ---------------------------------------------------------------

futile.logger::flog.info('\t- SortMeRNA')

sortMeRNA <- list()
sortMeRNA$files <- base::list.files(path = opt$input, full.names = T, recursive = T, pattern = '_rRNA_report.txt$')
sortMeRNA$data <- do.call(base::rbind, pbapply::pblapply(sortMeRNA$files, readSortMeRNA, cl = opt$threads))


# Preseq - Complexity curves ----------------------------------------------

futile.logger::flog.info('\t- PreSeq Complexity Curves')

preseq <- list()
preseq$files <- base::list.files(path = opt$input, full.names = T, recursive = T, pattern = 'ccurve.txt$')
preseq$data <- do.call(base::rbind, pbapply::pblapply(preseq$files, readPreseq, cl = opt$threads))


# Flagstat ----------------------------------------------------------------

futile.logger::flog.info('\t- Flagstat (Sambamba)')

flagstat <- list()
flagstat$files <- base::list.files(path = opt$input, full.names = T, recursive = T, pattern = '.flagstat$')
flagstat$data <- do.call(base::rbind, pbapply::pblapply(flagstat$files, readFlagstat.sambamba, cl = opt$threads))


# TIN ---------------------------------------------------------------------

futile.logger::flog.info('\t- Transcript Integity Number (TIN)')

TIN <- list()
TIN$files <- base::list.files(path = opt$input, full.names = T, recursive = T, pattern = '.summary.txt$')
TIN$data <- do.call(base::rbind, pbapply::pblapply(TIN$files, readTIN, cl = opt$threads))


# Generation of plots -----------------------------------------------------

futile.logger::flog.info('Generating interactive figures and tables.')

themePlots <- theme(
    legend.position = 'bottom',
    text = element_text(size = 10, family='Helvetica'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dotted'),
    panel.grid.minor = element_line(colour = 'grey80', linetype = 'dotted'),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
)


plots <- list()

# Overview of FASTQ files.
plots$overviewFastq <-
    TrimGalore$data %>% dplyr::group_by(sample, Direction) %>% dplyr::summarise(totalFile = base::length(base::unique(file))) %>%
    ggplot(., aes(x = sample, y = totalFile, fill = Direction)) +
    geom_bar(stat = 'identity', position = 'dodge', width = .9) +
    scale_fill_manual(values = c('#F47267', '#5880AB')) +
    labs(x = 'Samples', y = 'Number of input fastq files (R1 / R2)') +
    themePlots

# Overview of trimming - Reads.
plots$TrimGalore <-
    TrimGalore$data %>% dplyr::filter(grepl('Reads|Length', type, ignore.case = T)) %>%
    dplyr::mutate(type = factor(type, levels = c('Total reads processed', 'Reads with adapters', 'Below length cutoff (20 bp)', 'Reads written (passing filters)'))) %>%
    ggplot(., aes(x = file, y = value, color = type)) +
    geom_point() +
    scale_y_log10() +
    scale_color_manual(values = hues::iwanthue(4), guide = guide_legend(title = 'Categories', title.position = 'top', title.hjust = 0.5, nrow = 2, keywidth = 0.75, keyheight = 0.75)) +
    labs(x = 'Samples', y = 'Number of reads') +
    facet_wrap(. ~ sample, drop = T, scales = 'free_x') +
    themePlots + theme(axis.text.x = element_blank())

# SortMeRNA
plots$sortMeRNA <- sortMeRNA$data %>%
    dplyr::filter(type == 'Total reads passing E-value threshold') %>%
    ggplot(., aes(x = reorder(sample, -value.relative), y = value.relative, label = value)) +
    scale_y_continuous(labels = scales::percent) +
    geom_bar(stat = 'identity', position = 'dodge', fill = '#5880AB', color = 'black') +
    labs(x = 'Samples', y = 'Perc. of reads pessing E-value threshold\n(Are of rRNA origin)') +
    themePlots

# PreSeq - Complexity curves.
plots$preseq <-
    preseq$data %>%
    ggplot(., aes(x = TOTAL_READS, y = EXPECTED_DISTINCT, ymin = LOWER_0.95CI, ymax = UPPER_0.95CI, fill = sample, color = sample, group = sample)) +
    geom_line() +
    # geom_ribbon(alpha = 0.1, linetype = 'dotted') +
    scale_color_manual(values = hues::iwanthue(length(unique(preseq$data$sample)))) +
    scale_fill_manual(values = hues::iwanthue(length(unique(preseq$data$sample)))) +
    labs(x = 'Total Sequenced Reads', y = 'Expected Distinct Reads') +
    themePlots

# Flagstat.
plots$flagstat <- ggplot(flagstat$data, aes(x = reorder(type, -value), y = value, color = sample, group = sample)) +
    geom_boxplot(aes(group = type), color = 'grey75', alpha = .3, outlier.shape = NA) +
    geom_point() +
    geom_line(linetype = 'dotted', lwd = .1) +
    scale_y_log10() +
    scale_color_manual(values = hues::iwanthue(length(unique(flagstat$data$sample)))) +
    guides(color = F, group = F) +
    labs(x = 'Flagstat information', y = 'Number of reads') +
    themePlots

# TIN
plots$TIN <- ggplot(TIN$data, aes(x = reorder(sample, -`TIN(mean)`), y = `TIN(mean)`, group = sample)) +
    geom_bar(stat = 'identity', position = 'dodge', fill = '#5880AB', color = 'black') +
    labs(x = 'Samples', y = 'TIN (mean)') +
    themePlots
