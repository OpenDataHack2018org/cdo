/* Automatically created with makedoc, don't edit! */

std::vector<std::string> InfoHelp = {
  "NAME",
  "    info, infon, map - Information and simple statistics",
  "",
  "SYNOPSIS",
  "    <operator>  infiles",
  "",
  "DESCRIPTION",
  "    This module writes information about the structure and contents ",
  "    of all input files to standard output.  All input files need to have ",
  "    the same structure with the same variables on different timesteps.",
  "    The information displayed depends on the chosen operator.",
  "",
  "OPERATORS",
  "    info   Dataset information listed by parameter identifier",
  "           Prints information and simple statistics for each field of all",
  "           input datasets. For each field the operator prints one line "
  "with ",
  "           the following elements:",
  "           - Date and Time",
  "           - Level, Gridsize and number of Missing values",
  "           - Minimum, Mean and Maximum \\",
  "           The mean value is computed without the use of area weights!",
  "           - Parameter identifier",
  "    infon  Dataset information listed by parameter name",
  "           The same as operator info but using the name instead of the",
  "           identifier to label the parameter.",
  "    map    Dataset information and simple map",
  "           Prints information, simple statistics and a map for each field "
  "of all input",
  "           datasets. The map will be printed only for fields on a regular "
  "lon/lat grid.",
};

std::vector<std::string> SinfoHelp = {
  "NAME",
  "    sinfo, sinfon - Short information",
  "",
  "SYNOPSIS",
  "    <operator>  infiles",
  "",
  "DESCRIPTION",
  "    This module writes information about the structure of infiles to "
  "standard output.",
  "    infiles is an arbitrary number of input files. All input files need to "
  "have ",
  "    the same structure with the same variables on different timesteps.",
  "    The information displayed depends on the chosen operator.",
  "",
  "OPERATORS",
  "    sinfo   Short information listed by parameter identifier",
  "            Prints short information of a dataset. The information is "
  "divided into",
  "            4 sections. Section 1 prints one line per parameter with the "
  "following ",
  "            information:",
  "            - institute and source",
  "            - time c=constant v=varying",
  "            - type of statistical processing",
  "            - number of levels and z-axis number",
  "            - horizontal grid size and number",
  "            - data type",
  "            - parameter identifier",
  "            Section 2 and 3 gives a short overview of all grid and vertical "
  "coordinates.",
  "            And the last section contains short information of the time "
  "coordinate.",
  "    sinfon  Short information listed by parameter name",
  "            The same as operator sinfo but using the name instead of the ",
  "            identifier to label the parameter.",
};

std::vector<std::string> DiffHelp = {
  "NAME",
  "    diff, diffn - Compare two datasets field by field",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2",
  "",
  "DESCRIPTION",
  "    Compares the contents of two datasets field by field. The input "
  "datasets need to have the",
  "    same structure and its fields need to have the same header information "
  "and dimensions.",
  "",
  "OPERATORS",
  "    diff   Compare two datasets listed by parameter id",
  "           Provides statistics on differences between two datasets.",
  "           For each pair of fields the operator prints one line with the "
  "following information:",
  "           - Date and Time",
  "           - Level, Gridsize and number of Missing values",
  "           - Number of different values",
  "           - Occurrence of coefficient pairs with different signs (S)",
  "           - Occurrence of zero values (Z)",
  "           - Maxima of absolute difference of coefficient pairs",
  "           - Maxima of relative difference of non-zero coefficient pairs "
  "with equal signs",
  "           - Parameter identifier",
  "    diffn  Compare two datasets listed by parameter name",
  "           The same as operator diff. Using the name instead of the",
  "           identifier to label the parameter.",
};

std::vector<std::string> NinfoHelp = {
  "NAME",
  "    npar, nlevel, nyear, nmon, ndate, ntime, ngridpoints, ngrids - ",
  "    Print the number of parameters, levels or times",
  "",
  "SYNOPSIS",
  "    <operator>  infile",
  "",
  "DESCRIPTION",
  "    This module prints the number of variables, levels or times of the ",
  "    input dataset.",
  "",
  "OPERATORS",
  "    npar         Number of parameters",
  "                 Prints the number of parameters (variables).",
  "    nlevel       Number of levels",
  "                 Prints the number of levels for each variable.",
  "    nyear        Number of years",
  "                 Prints the number of different years.",
  "    nmon         Number of months",
  "                 Prints the number of different combinations of years and "
  "months.",
  "    ndate        Number of dates",
  "                 Prints the number of different dates.",
  "    ntime        Number of timesteps",
  "                 Prints the number of timesteps.",
  "    ngridpoints  Number of gridpoints",
  "                 Prints the number of gridpoints for each variable.",
  "    ngrids       Number of horizontal grids",
  "                 Prints the number of horizontal grids.",
};

std::vector<std::string> ShowinfoHelp = {
  "NAME",
  "    showformat, showcode, showname, showstdname, showatts, showattsglob, ",
  "    showlevel, showltype, showyear, showmon, showdate, showtime, "
  "showtimestamp - ",
  "    Show variables, levels or times",
  "",
  "SYNOPSIS",
  "    <operator>  infile",
  "",
  "DESCRIPTION",
  "    This module prints the format, variables, levels or times of the input "
  "dataset.",
  "",
  "OPERATORS",
  "    showformat     Show file format",
  "                   Prints the file format of the input dataset.",
  "    showcode       Show code numbers",
  "                   Prints the code number of all variables.",
  "    showname       Show variable names",
  "                   Prints the name of all variables.",
  "    showstdname    Show standard names",
  "                   Prints the standard name of all variables.",
  "    showatts       Show all attributes",
  "                   Prints all variable and global attributes.",
  "    showattsglob   Show all global attributes",
  "                   Prints all global attributes.",
  "    showlevel      Show levels",
  "                   Prints all levels for each variable.",
  "    showltype      Show GRIB level types",
  "                   Prints the GRIB level type for all z-axes.",
  "    showyear       Show years",
  "                   Prints all years.",
  "    showmon        Show months",
  "                   Prints all months.",
  "    showdate       Show date information",
  "                   Prints date information of all timesteps (format "
  "YYYY-MM-DD).",
  "    showtime       Show time information",
  "                   Prints time information of all timesteps (format "
  "hh:mm:ss).",
  "    showtimestamp  Show timestamp",
  "                   Prints timestamp of all timesteps (format "
  "YYYY-MM-DDThh:mm:ss).",
};

std::vector<std::string> ShowattributeHelp = {
  "NAME",
  "    showattribute, showattsvar - ",
  "    Show a global attribute, a variable attribute or all attributes of one "
  "variable",
  "",
  "SYNOPSIS",
  "    showattribute,attribute  infile",
  "    showattsvar[,var_nm]  infile",
  "",
  "DESCRIPTION",
  "    This operator prints attributes of a dataset.",
  "    If a global attribute should be printed, the attribute name can be "
  "specified as a parameter directly.",
  "    If a variable attribute should be printed, the following format is "
  "requested:",
  "    ",
  "      var_nm@att_nm",
  "    ",
  "       var_nm  Variable name. Example: pressure",
  "       att_nm  Attribute name. Example: units",
  "    ",
  "",
  "OPERATORS",
  "    showattribute  Show a global attribute or a variable attribute",
  "    showattsvar    Show all variable attributes.",
  "                   If var_nm is specified, only for a subset of variables.",
  "",
  "PARAMETER",
  "    attribute  STRING  Attribute in the format [var_nm@]att_nm",
  "    var_nm     STRING  Variable name",
};

std::vector<std::string> FiledesHelp = {
  "NAME",
  "    partab, codetab, griddes, zaxisdes, vct - Dataset description",
  "",
  "SYNOPSIS",
  "    <operator>  infile",
  "",
  "DESCRIPTION",
  "    This module provides operators to print meta information about a "
  "dataset.",
  "    The printed meta-data depends on the chosen operator.",
  "",
  "OPERATORS",
  "    partab    Parameter table",
  "              Prints all available meta information of the variables.",
  "    codetab   Parameter code table",
  "              Prints a code table with a description of all variables.",
  "              For each variable the operator prints one line listing the",
  "              code, name, description and units.",
  "    griddes   Grid description",
  "              Prints the description of all grids.",
  "    zaxisdes  Z-axis description",
  "              Prints the description of all z-axes.",
  "    vct       Vertical coordinate table",
  "              Prints the vertical coordinate table.",
};

std::vector<std::string> CopyHelp = {
  "NAME",
  "    copy, cat - Copy datasets",
  "",
  "SYNOPSIS",
  "    <operator>  infiles outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators to copy or concatenate datasets.",
  "    infiles is an arbitrary number of input files. All input files need to "
  "have ",
  "    the same structure with the same variables on different timesteps.",
  "",
  "OPERATORS",
  "    copy  Copy datasets",
  "          Copies all input datasets to outfile. ",
  "    cat   Concatenate datasets",
  "          Concatenates all input datasets and appends the result to the "
  "end ",
  "          of outfile. If outfile does not exist it will be created.",
};

std::vector<std::string> TeeHelp = {
  "NAME",
  "    tee - Duplicate a data stream",
  "",
  "SYNOPSIS",
  "    tee  infile outfile1 outfile2",
  "",
  "DESCRIPTION",
  "    This operator copies the input datasets to outfile1 and outfile2.",
  "    It can be used to store intermediate results to a file.",
};

std::vector<std::string> ReplaceHelp = {
  "NAME",
  "    replace - Replace variables",
  "",
  "SYNOPSIS",
  "    replace  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    The replace operator replaces variables in infile1 by variables from "
  "infile2 and write",
  "    the result to outfile. Both input datasets need to have the same number "
  "of timesteps.",
};

std::vector<std::string> DuplicateHelp = {
  "NAME",
  "    duplicate - Duplicates a dataset",
  "",
  "SYNOPSIS",
  "    duplicate[,ndup]  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator duplicates the contents of infile and writes the result "
  "to outfile.",
  "    The optional parameter sets the number of duplicates, the default is 2.",
  "",
  "PARAMETER",
  "    ndup  INTEGER  Number of duplicates, default is 2.",
};

std::vector<std::string> MergegridHelp = {
  "NAME",
  "    mergegrid - Merge grid",
  "",
  "SYNOPSIS",
  "    mergegrid  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Merges grid points of all variables from infile2 to infile1 and write "
  "the result to outfile.",
  "    Only the non missing values of infile2 will be used. The horizontal "
  "grid of infile2 should ",
  "    be smaller or equal to the grid of infile1 and the resolution must be "
  "the same.",
  "    Only rectilinear grids are supported. Both input files need to have the "
  "same variables ",
  "    and the same number of timesteps.",
};

std::vector<std::string> MergeHelp = {
  "NAME",
  "    merge, mergetime - Merge datasets",
  "",
  "SYNOPSIS",
  "    <operator>  infiles outfile",
  "",
  "DESCRIPTION",
  "    This module reads datasets from several input files, merges them and "
  "writes the resulting dataset to outfile.",
  "",
  "OPERATORS",
  "    merge      Merge datasets with different fields",
  "               Merges time series of different fields from several input "
  "datasets. The number ",
  "               of fields per timestep written to outfile is the sum of the "
  "field numbers ",
  "               per timestep in all input datasets. The time series on all "
  "input datasets are ",
  "               required to have different fields and the same number of "
  "timesteps.",
  "               The fields in each different input file either have to be "
  "different variables",
  "               or different levels of the same variable. A mixture of "
  "different variables on",
  "               different levels in different input files is not allowed.",
  "    mergetime  Merge datasets sorted by date and time",
  "               Merges all timesteps of all input files sorted by date and "
  "time.",
  "               All input files need to have the same structure with the "
  "same variables on ",
  "               different timesteps. After this operation every input "
  "timestep is in outfile ",
  "               and all timesteps are sorted by date and time.",
  "",
  "ENVIRONMENT",
  "    SKIP_SAME_TIME",
  "        If set to 1, skips all consecutive timesteps with a double entry of "
  "the same timestamp.",
  "",
  "NOTE",
  "    The operators in this module need to open all input files "
  "simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> SplitHelp = {
  "NAME",
  "    splitcode, splitparam, splitname, splitlevel, splitgrid, splitzaxis, ",
  "    splittabnum - Split a dataset",
  "",
  "SYNOPSIS",
  "    <operator>[,params]  infile obase",
  "",
  "DESCRIPTION",
  "    This module splits infile into pieces. The output files will be named "
  "<obase><xxx><suffix>",
  "    where suffix is the filename extension derived from the file format. "
  "xxx and the contents ",
  "    of the output files depends on the chosen operator. ",
  "    params is a comma separated list of processing parameters.",
  "",
  "OPERATORS",
  "    splitcode    Split code numbers",
  "                 Splits a dataset into pieces, one for each different code "
  "number.",
  "                 xxx will have three digits with the code number.",
  "    splitparam   Split parameter identifiers",
  "                 Splits a dataset into pieces, one for each different "
  "parameter identifier.",
  "                 xxx will be a string with the parameter identifier.",
  "    splitname    Split variable names",
  "                 Splits a dataset into pieces, one for each variable name.",
  "                 xxx will be a string with the variable name.",
  "    splitlevel   Split levels",
  "                 Splits a dataset into pieces, one for each different "
  "level.",
  "                 xxx will have six digits with the level.",
  "    splitgrid    Split grids",
  "                 Splits a dataset into pieces, one for each different grid.",
  "                 xxx will have two digits with the grid number.",
  "    splitzaxis   Split z-axes",
  "                 Splits a dataset into pieces, one for each different "
  "z-axis.",
  "                 xxx will have two digits with the z-axis number.",
  "    splittabnum  Split parameter table numbers",
  "                 Splits a dataset into pieces, one for each GRIB1 parameter "
  "table number.",
  "                 xxx will have three digits with the GRIB1 parameter table "
  "number.",
  "",
  "PARAMETER",
  "    swap            STRING  Swap the position of obase and xxx in the "
  "output filename",
  "    uuid=<attname>  STRING  Add a UUID as global attribute <attname> to "
  "each output file",
  "",
  "ENVIRONMENT",
  "    CDO_FILE_SUFFIX",
  "        Set the default file suffix. This suffix will be added to the "
  "output file ",
  "        names instead of the filename extension derived from the file "
  "format. ",
  "        Set this variable to NULL to disable the adding of a file suffix.",
  "",
  "NOTE",
  "    The operators in this module need to open all output files "
  "simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> SplittimeHelp = {
  "NAME",
  "    splithour, splitday, splitseas, splityear, splityearmon, splitmon - ",
  "    Split timesteps of a dataset",
  "",
  "SYNOPSIS",
  "    <operator>  infile obase",
  "    splitmon[,format]  infile obase",
  "",
  "DESCRIPTION",
  "    This module splits infile into  timesteps pieces. The output files will "
  "be named",
  "    <obase><xxx><suffix> where suffix is the filename extension derived "
  "from the file format. ",
  "    xxx and the contents of the output files depends on the chosen "
  "operator. ",
  "",
  "OPERATORS",
  "    splithour     Split hours",
  "                  Splits a file into pieces, one for each different hour.",
  "                  xxx will have two digits with the hour.",
  "    splitday      Split days",
  "                  Splits a file into pieces, one for each different day.",
  "                  xxx will have two digits with the day.",
  "    splitseas     Split seasons",
  "                  Splits a file into pieces, one for each different season.",
  "                  xxx will have three characters with the season.",
  "    splityear     Split years",
  "                  Splits a file into pieces, one for each different year.",
  "                  xxx will have four digits with the year (YYYY).",
  "    splityearmon  Split in years and months",
  "                  Splits a file into pieces, one for each different year "
  "and month.",
  "                  xxx will have six digits with the year and month "
  "(YYYYMM).",
  "    splitmon      Split months",
  "                  Splits a file into pieces, one for each different month.",
  "                  xxx will have two digits with the month.",
  "",
  "PARAMETER",
  "    format  STRING  C-style format for strftime() (e.g. %B for the full "
  "month name)",
  "",
  "ENVIRONMENT",
  "    CDO_FILE_SUFFIX",
  "        Set the default file suffix. This suffix will be added to the "
  "output file ",
  "        names instead of the filename extension derived from the file "
  "format. ",
  "        Set this variable to NULL to disable the adding of a file suffix.",
  "",
  "NOTE",
  "    The operators in this module need to open all output files "
  "simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> SplitselHelp = {
  "NAME",
  "    splitsel - Split selected timesteps",
  "",
  "SYNOPSIS",
  "    splitsel,nsets[,noffset[,nskip]]  infile obase",
  "",
  "DESCRIPTION",
  "    This operator splits infile into pieces, one for each adjacent",
  "    sequence t_1, ...., t_n of timesteps of the same selected time range.",
  "    The output files will be named <obase><nnnnnn><suffix> where nnnnnn is "
  "the ",
  "    sequence number and suffix is the filename extension derived from the "
  "file format.",
  "",
  "PARAMETER",
  "    nsets    INTEGER  Number of input timesteps for each output file",
  "    noffset  INTEGER  Number of input timesteps skipped before the first "
  "timestep range (optional)",
  "    nskip    INTEGER  Number of input timesteps skipped between timestep "
  "ranges (optional)",
  "",
  "ENVIRONMENT",
  "    CDO_FILE_SUFFIX",
  "        Set the default file suffix. This suffix will be added to the "
  "output file ",
  "        names instead of the filename extension derived from the file "
  "format. ",
  "        Set this variable to NULL to disable the adding of a file suffix.",
};

std::vector<std::string> DistgridHelp = {
  "NAME",
  "    distgrid - Distribute horizontal grid",
  "",
  "SYNOPSIS",
  "    distgrid,nx[,ny]  infile obase",
  "",
  "DESCRIPTION",
  "    This operator distributes a dataset into smaller pieces. Each output "
  "file contains a different region of the",
  "    horizontal source grid. A target grid region contains a structured "
  "longitude/latitude box of the source grid.",
  "    Only rectilinear and curvilinear source grids are supported by this "
  "operator.",
  "    The number of different regions can be specified with the parameter nx "
  "and ny. The output files will be named ",
  "    <obase><xxx><suffix> where suffix is the filename extension derived "
  "from the file format. xxx will have five digits with ",
  "    the number of the target region.",
  "",
  "PARAMETER",
  "    nx  INTEGER  Number of regions in x direction",
  "    ny  INTEGER  Number of regions in y direction [default: 1]",
  "",
  "NOTE",
  "    This operator needs to open all output files simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> CollgridHelp = {
  "NAME",
  "    collgrid - Collect horizontal grid",
  "",
  "SYNOPSIS",
  "    collgrid[,nx[,names]]  infiles outfile",
  "",
  "DESCRIPTION",
  "    This operator collects the data of the input files to one output file. ",
  "    All input files need to have the same variables and the same number of "
  "timesteps on a different",
  "    horizonal grid region. A source region must be a structured "
  "longitude/latitude grid box.",
  "    The parameter nx needs to be specified only for non regular lon/lat "
  "grids.",
  "",
  "PARAMETER",
  "    nx     INTEGER  Number of regions in x direction [default: number of "
  "input files]",
  "    names  STRING   Comma separated list of variable names [default: all "
  "variables]",
  "",
  "NOTE",
  "    This operator needs to open all input files simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> SelectHelp = {
  "NAME",
  "    select, delete - Select fields",
  "",
  "SYNOPSIS",
  "    <operator>,params  infiles outfile",
  "",
  "DESCRIPTION",
  "    This module selects some fields from infiles and writes them to "
  "outfile.",
  "    infiles is an arbitrary number of input files. All input files need to "
  "have ",
  "    the same structure with the same variables on different timesteps.",
  "    The fields selected depends on the chosen parameters. Parameter is a "
  "comma",
  "    separated list of \"key=value\" pairs. Wildcards can be used for string "
  "parameter.",
  "",
  "OPERATORS",
  "    select  Select fields",
  "            Selects all fields with parameters in a user given list.",
  "    delete  Delete fields",
  "            Deletes all fields with parameters in a user given list.",
  "",
  "PARAMETER",
  "    name              STRING  Comma separated list of variable names.",
  "    param             STRING  Comma separated list of parameter "
  "identifiers.",
  "    code              INTEGER Comma separated list of code numbers.",
  "    level             FLOAT   Comma separated list of vertical levels.",
  "    levidx            INTEGER Comma separated list of index of levels.",
  "    zaxisname         STRING  Comma separated list of zaxis names.",
  "    zaxisnum          INTEGER Comma separated list of zaxis numbers.",
  "    ltype             INTEGER Comma separated list of GRIB level types.",
  "    gridname          STRING  Comma separated list of grid names.",
  "    gridnum           INTEGER Comma separated list of grid numbers.",
  "    steptype          STRING  Comma separated list of timestep types.",
  "    date              STRING  Comma separated list of dates (format "
  "YYYY-MM-DDThh:mm:ss).",
  "    startdate         STRING  Start date (format YYYY-MM-DDThh:mm:ss).",
  "    enddate           STRING  End date (format YYYY-MM-DDThh:mm:ss).",
  "    minute            INTEGER Comma separated list of minutes.",
  "    hour              INTEGER Comma separated list of hours.",
  "    day               INTEGER Comma separated list of days.",
  "    month             INTEGER Comma separated list of months.",
  "    season            STRING  Comma separated list of seasons (substring of "
  "DJFMAMJJASOND or ANN).",
  "    year              INTEGER Comma separated list of years.",
  "    timestep          INTEGER Comma separated list of timesteps. Negative "
  "values selects timesteps from the end (NetCDF only).",
  "    timestep_of_year  INTEGER Comma separated list of timesteps of year.",
  "    timestepmask      STRING  Read timesteps from a mask file.",
};

std::vector<std::string> SelmultiHelp = {
  "NAME",
  "    selmulti, delmulti, changemulti - Select multiple fields via GRIB1 "
  "parameters",
  "",
  "SYNOPSIS",
  "    <operator>,selection-specification  infile outfile",
  "",
  "DESCRIPTION",
  "    This module selects multiple fields from infile and writes them to "
  "outfile.",
  "    selection-specification is a filename or in-place string with the "
  "selection specification.",
  "    Each selection-specification has the following compact notation "
  "format: ",
  "    ",
  "       <type>(parameters; leveltype(s); levels)",
  "    ",
  "    type      "
  "    sel for select or del for delete (optional)",
  "    parameters"
  "    GRIB1 parameter code number",
  "    leveltype "
  "    GRIB1 level type",
  "    levels    "
  "    value of each level",
  "    ",
  "    Examples:",
  "    ",
  "       (1; 103; 0) ",
  "       (33,34; 105; 10) ",
  "       (11,17; 105; 2) ",
  "       (71,73,74,75,61,62,65,117,67,122,121,11,131,66,84,111,112; 105; 0) ",
  "    ",
  "    The following descriptive notation can also be used for selection "
  "specification from a file:",
  "    ",
  "       SELECT/DELETE, PARAMETER=parameters, LEVTYPE=leveltye(s), "
  "LEVEL=levels",
  "    ",
  "    Examples:",
  "    ",
  "       SELECT, PARAMETER=1, LEVTYPE=103, LEVEL=0 ",
  "       SELECT, PARAMETER=33/34, LEVTYPE=105, LEVEL=10 ",
  "       SELECT, PARAMETER=11/17, LEVTYPE=105, LEVEL=2 ",
  "       SELECT, PARAMETER=71/73/74/75/61/62/65/117/67/122, LEVTYPE=105, "
  "LEVEL=0 ",
  "       DELETE, PARAMETER=128, LEVTYPE=109, LEVEL=* ",
  "    ",
  "    The following will convert Pressure from Pa into hPa; Temp from Kelvin "
  "to Celsius: ",
  "       SELECT, PARAMETER=1, LEVTYPE= 103, LEVEL=0, SCALE=0.01 ",
  "       SELECT, PARAMETER=11, LEVTYPE=105, LEVEL=2, OFFSET=273.15 ",
  "    If SCALE and/or OFFSET are defined, then the data values are scaled as "
  "SCALE*(VALUE-OFFSET).",
  "",
  "OPERATORS",
  "    selmulti     Select multiple fields",
  "    delmulti     Delete multiple fields",
  "    changemulti  Change identication of multiple fields",
};

std::vector<std::string> SelvarHelp = {
  "NAME",
  "    selparam, delparam, selcode, delcode, selname, delname, selstdname, "
  "sellevel, ",
  "    sellevidx, selgrid, selzaxis, selzaxisname, selltype, seltabnum - "
  "Select fields",
  "",
  "SYNOPSIS",
  "    <operator>,params  infile outfile",
  "    selcode,codes  infile outfile",
  "    delcode,codes  infile outfile",
  "    selname,names  infile outfile",
  "    delname,names  infile outfile",
  "    selstdname,stdnames  infile outfile",
  "    sellevel,levels  infile outfile",
  "    sellevidx,levidx  infile outfile",
  "    selgrid,grids  infile outfile",
  "    selzaxis,zaxes  infile outfile",
  "    selzaxisname,zaxisnames  infile outfile",
  "    selltype,ltypes  infile outfile",
  "    seltabnum,tabnums  infile outfile",
  "",
  "DESCRIPTION",
  "    This module selects some fields from infile and writes them to outfile.",
  "    The fields selected depends on the chosen operator and the parameters.",
  "",
  "OPERATORS",
  "    selparam      Select parameters by identifier",
  "                  Selects all fields with parameter identifiers in a user "
  "given list.",
  "    delparam      Delete parameters by identifier",
  "                  Deletes all fields with parameter identifiers in a user "
  "given list.",
  "    selcode       Select parameters by code number",
  "                  Selects all fields with code numbers in a user given "
  "list.",
  "    delcode       Delete parameters by code number",
  "                  Deletes all fields with code numbers in a user given "
  "list.",
  "    selname       Select parameters by name",
  "                  Selects all fields with parameter names in a user given "
  "list.",
  "    delname       Delete parameters by name",
  "                  Deletes all fields with parameter names in a user given "
  "list.",
  "    selstdname    Select parameters by standard name",
  "                  Selects all fields with standard names in a user given "
  "list.",
  "    sellevel      Select levels",
  "                  Selects all fields with levels in a user given list.",
  "    sellevidx     Select levels by index",
  "                  Selects all fields with index of levels in a user given "
  "list.",
  "    selgrid       Select grids",
  "                  Selects all fields with grids in a user given list.",
  "    selzaxis      Select z-axes",
  "                  Selects all fields with z-axes in a user given list.",
  "    selzaxisname  Select z-axes by name",
  "                  Selects all fields with z-axis names in a user given "
  "list.",
  "    selltype      Select GRIB level types",
  "                  Selects all fields with GRIB level type in a user given "
  "list.",
  "    seltabnum     Select parameter table numbers",
  "                  Selects all fields with parameter table numbers in a user "
  "given list.",
  "",
  "PARAMETER",
  "    params      INTEGER  Comma separated list of parameter identifiers",
  "    codes       INTEGER  Comma separated list of code numbers",
  "    names       STRING   Comma separated list of variable names",
  "    stdnames    STRING   Comma separated list of standard names",
  "    levels      FLOAT    Comma separated list of vertical levels",
  "    levidx      INTEGER  Comma separated list of index of levels",
  "    ltypes      INTEGER  Comma separated list of GRIB level types",
  "    grids       STRING   Comma separated list of grid names or numbers",
  "    zaxes       STRING   Comma separated list of z-axis types or numbers",
  "    zaxisnames  STRING   Comma separated list of z-axis names",
  "    tabnums     INTEGER  Comma separated list of parameter table numbers",
};

std::vector<std::string> SeltimeHelp = {
  "NAME",
  "    seltimestep, seltime, selhour, selday, selmonth, selyear, selseason, "
  "seldate, ",
  "    selsmon - Select timesteps",
  "",
  "SYNOPSIS",
  "    seltimestep,timesteps  infile outfile",
  "    seltime,times  infile outfile",
  "    selhour,hours  infile outfile",
  "    selday,days  infile outfile",
  "    selmonth,months  infile outfile",
  "    selyear,years  infile outfile",
  "    selseason,seasons  infile outfile",
  "    seldate,startdate[,enddate]  infile outfile",
  "    selsmon,month[,nts1[,nts2]]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module selects user specified timesteps from infile and writes "
  "them to outfile.",
  "    The timesteps selected depends on the chosen operator and the "
  "parameters.",
  "",
  "OPERATORS",
  "    seltimestep  Select timesteps",
  "                 Selects all timesteps with a timestep in a user given "
  "list.",
  "    seltime      Select times",
  "                 Selects all timesteps with a time in a user given list.",
  "    selhour      Select hours",
  "                 Selects all timesteps with a hour in a user given list.",
  "    selday       Select days",
  "                 Selects all timesteps with a day in a user given list.",
  "    selmonth     Select months",
  "                 Selects all timesteps with a month in a user given list.",
  "    selyear      Select years",
  "                 Selects all timesteps with a year in a user given list.",
  "    selseason    Select seasons",
  "                 Selects all timesteps with a month of a season in a user "
  "given list.",
  "    seldate      Select dates",
  "                 Selects all timesteps with a date in a user given range.",
  "    selsmon      Select single month",
  "                 Selects a month and optional an arbitrary number of "
  "timesteps before and after this month.",
  "",
  "PARAMETER",
  "    timesteps  INTEGER  Comma separated list of timesteps. Negative values "
  "selects timesteps from the end (NetCDF only).",
  "    times      STRING   Comma separated list of times (format hh:mm:ss).",
  "    hours      INTEGER  Comma separated list of hours.",
  "    days       INTEGER  Comma separated list of days.",
  "    months     INTEGER  Comma separated list of months.",
  "    years      INTEGER  Comma separated list of years.",
  "    seasons    STRING   Comma separated list of seasons (substring of "
  "DJFMAMJJASOND or ANN).",
  "    startdate  STRING   Start date (format YYYY-MM-DDThh:mm:ss).",
  "    enddate    STRING   End date (format YYYY-MM-DDThh:mm:ss) [default: "
  "startdate].",
  "    nts1       INTEGER  Number of timesteps before the selected month "
  "[default: 0].",
  "    nts2       INTEGER  Number of timesteps after the selected month "
  "[default: nts1].",
};

std::vector<std::string> SelboxHelp = {
  "NAME",
  "    sellonlatbox, selindexbox - Select a box of a field",
  "",
  "SYNOPSIS",
  "    sellonlatbox,lon1,lon2,lat1,lat2  infile outfile",
  "    selindexbox,idx1,idx2,idy1,idy2  infile outfile",
  "",
  "DESCRIPTION",
  "    Selects a box of the rectangularly understood field.",
  "",
  "OPERATORS",
  "    sellonlatbox  Select a longitude/latitude box",
  "                  Selects a regular longitude/latitude box. The user has to "
  "give the longitudes and latitudes of the ",
  "                  edges of the box. Considered are only those grid cells "
  "with the grid center inside the lon/lat box.",
  "                  For rotated lon/lat grids the parameter needs to be "
  "rotated coordinates.",
  "    selindexbox   Select an index box",
  "                  Selects an index box. The user has to give the indexes of "
  "the edges of the box. The index of the ",
  "                  left edge may be greater then that of the right edge.",
  "",
  "PARAMETER",
  "    lon1  FLOAT    Western longitude",
  "    lon2  FLOAT    Eastern longitude",
  "    lat1  FLOAT    Southern or northern latitude",
  "    lat2  FLOAT    Northern or southern latitude",
  "    idx1  INTEGER  Index of first longitude (1 - nlon)",
  "    idx2  INTEGER  Index of last longitude (1 - nlon)",
  "    idy1  INTEGER  Index of first latitude (1 - nlat)",
  "    idy2  INTEGER  Index of last latitude (1 - nlat)",
};

std::vector<std::string> SelgridcellHelp = {
  "NAME",
  "    selgridcell, delgridcell - Select grid cells",
  "",
  "SYNOPSIS",
  "    <operator>,indexes  infile outfile",
  "",
  "DESCRIPTION",
  "    Selects grid cells of all fields from infile. The user has to give the "
  "indexes of each",
  "    grid cell. The resulting grid in outfile is unstructured.",
  "",
  "OPERATORS",
  "    selgridcell  Select grid cells",
  "    delgridcell  Delete grid cells",
  "",
  "PARAMETER",
  "    indexes  INTEGER  Comma separated list of indexes",
};

std::vector<std::string> SamplegridHelp = {
  "NAME",
  "    samplegrid - Resample grid",
  "",
  "SYNOPSIS",
  "    samplegrid,factor  infile outfile",
  "",
  "DESCRIPTION",
  "    This is a special operator for resampling the horizontal grid.",
  "    No interpolation takes place. Resample factor=2 means every second grid "
  "point is removed.",
  "    Only rectilinear and curvilinear source grids are supported by this "
  "operator.",
  "",
  "PARAMETER",
  "    factor  INTEGER  Resample factor, typically 2, which will half the "
  "resolution",
};

std::vector<std::string> CondHelp = {
  "NAME",
  "    ifthen, ifnotthen - Conditional select one field",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module selects field elements from infile2 with respect to infile1 "
  "and writes them ",
  "    to outfile. The fields in infile1 are handled as a mask. A value ",
  "    not equal to zero is treated as \"true\", zero is treated as \"false\".",
  "    The number of fields in infile1 has either to be the same as in infile2 "
  "or the",
  "    same as in one timestep of infile2 or only one.",
  "    The fields in outfile inherit the meta data from infile2.",
  "",
  "OPERATORS",
  "    ifthen     If then",
  "                        / i_2(t,x) if i_1([t,]x) NE 0  AND  i_1([t,]x) NE "
  "miss",
  "               o(t,x) =",
  "                        \\ miss     if i_1([t,]x) EQ 0  OR   i_1([t,]x) EQ "
  "miss",
  "    ifnotthen  If not then",
  "                        / i_2(t,x) if i_1([t,]x) EQ 0  AND  i_1([t,]x) NE "
  "miss",
  "               o(t,x) = ",
  "                        \\ miss     if i_1([t,]x) NE 0  OR   i_1([t,]x) EQ "
  "miss",
};

std::vector<std::string> Cond2Help = {
  "NAME",
  "    ifthenelse - Conditional select  two fields",
  "",
  "SYNOPSIS",
  "    ifthenelse  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator selects field elements from infile2 or infile3 with "
  "respect to",
  "    infile1 and writes them to outfile. The fields in infile1 are handled "
  "as a mask.",
  "    A value not equal to zero is treated as \"true\", zero is treated as "
  "\"false\".",
  "    The number of fields in infile1 has either to be the same as in infile2 "
  "or the ",
  "    same as in one timestep of infile2 or only one.",
  "    infile2 and infile3 need to have the same number of fields.",
  "    The fields in outfile inherit the meta data from infile2.",
  "    ",
  "              / i_2(t,x) if i_1([t,]x) NE 0  AND  i_1([t,]x) NE miss",
  "    o(t,x) = <  i_3(t,x) if i_1([t,]x) EQ 0  AND  i_1([t,]x) NE miss",
  "              \\ miss     if i_1([t,]x) EQ miss",
};

std::vector<std::string> CondcHelp = {
  "NAME",
  "    ifthenc, ifnotthenc - Conditional select a constant",
  "",
  "SYNOPSIS",
  "    <operator>,c  infile outfile",
  "",
  "DESCRIPTION",
  "    This module creates fields with a constant value or missing value.",
  "    The fields in infile are handled as a mask. A value not equal ",
  "    to zero is treated as \"true\", zero is treated as \"false\".",
  "",
  "OPERATORS",
  "    ifthenc     If then constant",
  "                         / c      if i(t,x) NE 0  AND  i(t,x) NE miss",
  "                o(t,x) =",
  "                         \\ miss   if i(t,x) EQ 0  OR   i(t,x) EQ miss",
  "    ifnotthenc  If not then constant",
  "                         / c      if i(t,x) EQ 0  AND  i(t,x) NE miss",
  "                o(t,x) =",
  "                         \\ miss   if i(t,x) NE 0  OR   i(t,x) EQ miss",
  "",
  "PARAMETER",
  "    c  FLOAT  Constant",
};

std::vector<std::string> MapReduceHelp = {
  "NAME",
  "    reducegrid - Reduce fields to user-defined mask",
  "",
  "SYNOPSIS",
  "    reducegrid,mask[,limitCoordsOutput]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module holds an operator for data reduction based on a user "
  "defined mask.",
  "    The output grid is unstructured and includes coordinate bounds. Bounds "
  "can be",
  "    avoided by using the additional 'nobounds' keyword. With 'nocoords' "
  "given,",
  "    coordinates a completely suppressed.",
  "",
  "PARAMETER",
  "    mask               STRING file which holds the mask field",
  "    limitCoordsOutput  STRING optional parameter to limit coordinates "
  "output: 'nobounds' disables coordinate bounds, 'nocoords' avoids all "
  "coordinate information",
};

std::vector<std::string> CompHelp = {
  "NAME",
  "    eq, ne, le, lt, ge, gt - Comparison of two fields",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module compares two datasets field by field. The resulting",
  "    field is a mask containing 1 if the comparison is true and 0 if not. ",
  "    The number of fields in infile1 should be the same as in infile2.",
  "    One of the input files can contain only one timestep or one field.",
  "    The fields in outfile inherit the meta data from infile1 or infile2.",
  "    The type of comparison depends on the chosen operator.",
  "",
  "OPERATORS",
  "    eq  Equal",
  "                  /   1   if i_1(t,x) EQ i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "        o(t,x) = <    0   if i_1(t,x) NE i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "                  \\  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss",
  "    ne  Not equal",
  "                  /   1   if i_1(t,x) NE i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "        o(t,x) = <    0   if i_1(t,x) EQ i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "                  \\  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss",
  "    le  Less equal",
  "                  /   1   if i_1(t,x) LE i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "        o(t,x) = <    0   if i_1(t,x) GT i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "                  \\  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss",
  "    lt  Less than",
  "                  /   1   if i_1(t,x) LT i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "        o(t,x) = <    0   if i_1(t,x) GE i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "                  \\  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss",
  "    ge  Greater equal",
  "                  /   1   if i_1(t,x) GE i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "        o(t,x) = <    0   if i_1(t,x) LT i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "                  \\  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss",
  "    gt  Greater than",
  "                  /   1   if i_1(t,x) GT i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "        o(t,x) = <    0   if i_1(t,x) LE i_2(t,x)  AND  i_1(t,x),i_2(t,x) "
  "NE miss",
  "                  \\  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss",
};

std::vector<std::string> CompcHelp = {
  "NAME",
  "    eqc, nec, lec, ltc, gec, gtc - Comparison of a field with a constant",
  "",
  "SYNOPSIS",
  "    <operator>,c  infile outfile",
  "",
  "DESCRIPTION",
  "    This module compares all fields of a dataset with a constant. The "
  "resulting",
  "    field is a mask containing 1 if the comparison is true and 0 if not.",
  "    The type of comparison depends on the chosen operator.",
  "",
  "OPERATORS",
  "    eqc  Equal constant",
  "                   /   1   if i(t,x) EQ c     AND  i(t,x),c NE miss",
  "         o(t,x) = <    0   if i(t,x) NE c     AND  i(t,x),c NE miss",
  "                   \\  miss if i(t,x) EQ miss  OR   c EQ miss",
  "    nec  Not equal constant",
  "                   /   1   if i(t,x) NE c     AND  i(t,x),c NE miss",
  "         o(t,x) = <    0   if i(t,x) EQ c     AND  i(t,x),c NE miss",
  "                   \\  miss if i(t,x) EQ miss  OR   c EQ miss",
  "    lec  Less equal constant",
  "                   /   1   if i(t,x) LE c     AND  i(t,x),c NE miss",
  "         o(t,x) = <    0   if i(t,x) GT c     AND  i(t,x),c NE miss",
  "                   \\  miss if i(t,x) EQ miss  OR   c EQ miss",
  "    ltc  Less than constant",
  "                   /   1   if i(t,x) LT c     AND  i(t,x),c NE miss",
  "         o(t,x) = <    0   if i(t,x) GE c     AND  i(t,x),c NE miss",
  "                   \\  miss if i(t,x) EQ miss  OR   c EQ miss",
  "    gec  Greater equal constant",
  "                   /   1   if i(t,x) GE c     AND  i(t,x),c NE miss",
  "         o(t,x) = <    0   if i(t,x) LT c     AND  i(t,x),c NE miss",
  "                   \\  miss if i(t,x) EQ miss  OR   c EQ miss",
  "    gtc  Greater than constant",
  "                   /   1   if i(t,x) GT c     AND  i(t,x),c NE miss",
  "         o(t,x) = <    0   if i(t,x) LE c     AND  i(t,x),c NE miss",
  "                   \\  miss if i(t,x) EQ miss  OR   c EQ miss",
  "",
  "PARAMETER",
  "    c  FLOAT  Constant",
};

std::vector<std::string> SetattributeHelp = {
  "NAME",
  "    setattribute - Set attributes",
  "",
  "SYNOPSIS",
  "    setattribute,attributes  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator sets attributes of a dataset. Each attribute has the "
  "following structure:",
  "    ",
  "      [var_nm@]att_nm=att_val",
  "    ",
  "       var_nm  Variable name (optional). Example: pressure",
  "       att_nm  Attribute name. Example: units",
  "       att_val Comma separated list of attribute values. Example: pascal",
  "    ",
  "    The value of var_nm is the name of the variable containing the "
  "attribute (named att_nm) that",
  "    you want to set. Use wildcards to set the attribute att_nm to more than "
  "one variable.",
  "    A value of var_nm of '*' will set the attribute att_nm to all data "
  "variables.",
  "    If var_nm is missing then att_nm refers to a global attribute.",
  "    ",
  "    The value of att_nm is the name of the attribute you want to set.",
  "    ",
  "    The value of att_val is the contents of the attribute att_nm. att_val "
  "may be a single value",
  "    or one-dimensional array of elements. The type of the attribute value "
  "will be detected",
  "    automaticly from the contents of the value.",
  "    ",
  "    A special meaning has the attribute name FILE. If this is the 1st "
  "attribute then all attributes",
  "    are read from a file specified in the value of att_val.",
  "",
  "PARAMETER",
  "    attributes  STRING  Comma separated list of attributes. ",
};

std::vector<std::string> SetpartabHelp = {
  "NAME",
  "    setpartabp, setpartabn - Set parameter table",
  "",
  "SYNOPSIS",
  "    <operator>,table[,convert]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module transforms data and metadata of infile via a parameter "
  "table and writes the result to outfile.",
  "    A parameter table is an ASCII formatted file with a set of parameter "
  "entries for each variable. Each new set have to",
  "    start with \"\\&parameter\" and to end with \"/\".",
  "    ",
  "    The following parameter table entries are supported:",
  "    ",
  "     Entry           & Type        & Description      ",
  "     name            & WORD        & Name of the variable",
  "     out_name        & WORD        & New name of the variable",
  "     param           & WORD        & Parameter identifier (GRIB1: "
  "code[.tabnum];  GRIB2: num[.cat[.dis]])",
  "     out_param       & WORD        & New parameter identifier",
  "     type            & WORD        & Data type (real or double)",
  "     standard_name   & WORD        & As defined in the CF standard name "
  "table",
  "     long_name       & STRING      & Describing the variable",
  "     units           & STRING      & Specifying the units for the variable",
  "     comment         & STRING      & Information concerning the variable",
  "     cell_methods    & STRING      & Information concerning calculation of "
  "means or climatologies",
  "     cell_measures   & STRING      & Indicates the names of the variables "
  "containing cell areas and volumes",
  "     missing_value   & FLOAT       & Specifying how missing data will be "
  "identified",
  "     valid_min       & FLOAT       & Minimum valid value",
  "     valid_max       & FLOAT       & Maximum valid value",
  "     ok_min_mean_abs & FLOAT       & Minimum absolute mean",
  "     ok_max_mean_abs & FLOAT       & Maximum absolute mean",
  "     factor          & FLOAT       & Scale factor",
  "     delete          & INTEGER     & Set to 1 to delete variable",
  "     convert         & INTEGER     & Set to 1 to convert the unit if "
  "necessary",
  "    ",
  "    Unsupported parameter table entries are stored as variable attributes.",
  "    The search key for the variable depends on the operator. Use setpartabn "
  "to search variables by the name.",
  "    This is typically used for NetCDF datasets. The operator setpartabp "
  "searches variables by the parameter ID.",
  "",
  "OPERATORS",
  "    setpartabp  Set parameter table",
  "                Search variables by the parameter identifier.",
  "    setpartabn  Set parameter table",
  "                Search variables by name.",
  "",
  "PARAMETER",
  "    table    STRING   Parameter table file or name",
  "    convert  STRING   Converts the units if necessary",
};

std::vector<std::string> SetHelp = {
  "NAME",
  "    setcodetab, setcode, setparam, setname, setunit, setlevel, setltype - ",
  "    Set field info",
  "",
  "SYNOPSIS",
  "    setcodetab,table  infile outfile",
  "    setcode,code  infile outfile",
  "    setparam,param  infile outfile",
  "    setname,name  infile outfile",
  "    setunit,unit  infile outfile",
  "    setlevel,level  infile outfile",
  "    setltype,ltype  infile outfile",
  "",
  "DESCRIPTION",
  "    This module sets some field information. Depending on the chosen "
  "operator the ",
  "    parameter table, code number, parameter identifier, variable name or "
  "level is set.",
  "",
  "OPERATORS",
  "    setcodetab  Set parameter code table",
  "                Sets the parameter code table for all variables.",
  "    setcode     Set code number",
  "                Sets the code number for all variables to the same given "
  "value.",
  "    setparam    Set parameter identifier",
  "                Sets the parameter identifier of the first variable.",
  "    setname     Set variable name",
  "                Sets the name of the first variable.",
  "    setunit     Set variable unit",
  "                Sets the unit of the first variable.",
  "    setlevel    Set level",
  "                Sets the first level of all variables.",
  "    setltype    Set GRIB level type",
  "                Sets the GRIB level type of all variables.",
  "",
  "PARAMETER",
  "    table  STRING   Parameter table file or name",
  "    code   INTEGER  Code number",
  "    param  STRING   Parameter identifier (GRIB1: code[.tabnum]; GRIB2: "
  "num[.cat[.dis]])",
  "    name   STRING   Variable name",
  "    level  FLOAT    New level",
  "    ltype  INTEGER  GRIB level type",
};

std::vector<std::string> SettimeHelp = {
  "NAME",
  "    setdate, settime, setday, setmon, setyear, settunits, settaxis, "
  "settbounds, ",
  "    setreftime, setcalendar, shifttime - Set time",
  "",
  "SYNOPSIS",
  "    setdate,date  infile outfile",
  "    settime,time  infile outfile",
  "    setday,day  infile outfile",
  "    setmon,month  infile outfile",
  "    setyear,year  infile outfile",
  "    settunits,units  infile outfile",
  "    settaxis,date,time[,inc]  infile outfile",
  "    settbounds,frequency  infile outfile",
  "    setreftime,date,time[,units]  infile outfile",
  "    setcalendar,calendar  infile outfile",
  "    shifttime,sval  infile outfile",
  "",
  "DESCRIPTION",
  "    This module sets the time axis or part of the time axis. Which part of "
  "the",
  "    time axis is overwritten/created depends on the chosen operator.",
  "",
  "OPERATORS",
  "    setdate      Set date",
  "                 Sets the date in every timestep to the same given value.",
  "    settime      Set time of the day",
  "                 Sets the time in every timestep to the same given value.",
  "    setday       Set day",
  "                 Sets the day in every timestep to the same given value.",
  "    setmon       Set month",
  "                 Sets the month in every timestep to the same given value.",
  "    setyear      Set year",
  "                 Sets the year in every timestep to the same given value.",
  "    settunits    Set time units",
  "                 Sets the base units of a relative time axis.",
  "    settaxis     Set time axis",
  "                 Sets the time axis.",
  "    settbounds   Set time bounds",
  "                 Sets the time bounds.",
  "    setreftime   Set reference time",
  "                 Sets the reference time of a relative time axis.",
  "    setcalendar  Set calendar",
  "                 Sets the calendar of a relative time axis.",
  "    shifttime    Shift timesteps",
  "                 Shifts all timesteps by the parameter sval.",
  "",
  "PARAMETER",
  "    day        INTEGER  Value of the new day",
  "    month      INTEGER  Value of the new month",
  "    year       INTEGER  Value of the new year",
  "    units      STRING   Base units of the time axis (seconds, minutes, "
  "hours, days, months, years)",
  "    date       STRING   Date (format: YYYY-MM-DD)",
  "    time       STRING   Time (format: hh:mm:ss)",
  "    inc        STRING   Optional increment (seconds, minutes, hours, days, "
  "months, years) [default: 1hour]",
  "    frequency  STRING   Frequency of the time series (hour, day, month, "
  "year)",
  "    calendar   STRING   Calendar (standard, proleptic_gregorian, 360_day, "
  "365_day, 366_day)",
  "    sval       STRING   Shift value (e.g. -3hour)",
};

std::vector<std::string> ChangeHelp = {
  "NAME",
  "    chcode, chparam, chname, chunit, chlevel, chlevelc, chlevelv - ",
  "    Change field header",
  "",
  "SYNOPSIS",
  "    chcode,oldcode,newcode[,...]  infile outfile",
  "    chparam,oldparam,newparam,...  infile outfile",
  "    chname,oldname,newname,...  infile outfile",
  "    chunit,oldunit,newunit,...  infile outfile",
  "    chlevel,oldlev,newlev,...  infile outfile",
  "    chlevelc,code,oldlev,newlev  infile outfile",
  "    chlevelv,name,oldlev,newlev  infile outfile",
  "",
  "DESCRIPTION",
  "    This module reads fields from infile, changes some header values",
  "    and writes the results to outfile. The kind of changes depends on ",
  "    the chosen operator.",
  "",
  "OPERATORS",
  "    chcode    Change code number",
  "              Changes some user given code numbers to new user given "
  "values.",
  "    chparam   Change parameter identifier",
  "              Changes some user given parameter identifiers to new user "
  "given values.",
  "    chname    Change variable name",
  "              Changes some user given variable names to new user given "
  "names.",
  "    chunit    Change variable unit",
  "              Changes some user given variable units to new user given "
  "units.",
  "    chlevel   Change level",
  "              Changes some user given levels to new user given values.",
  "    chlevelc  Change level of one code",
  "              Changes one level of a user given code number.",
  "    chlevelv  Change level of one variable",
  "              Changes one level of a user given variable name.",
  "",
  "PARAMETER",
  "    code                   INTEGER  Code number",
  "    oldcode,newcode,...    INTEGER  Pairs of old and new code numbers",
  "    oldparam,newparam,...  STRING   Pairs of old and new parameter "
  "identifiers",
  "    name                   STRING   Variable name",
  "    oldname,newname,...    STRING   Pairs of old and new variable names",
  "    oldlev                 FLOAT    Old level",
  "    newlev                 FLOAT    New level",
  "    oldlev,newlev,...      FLOAT    Pairs of old and new levels",
};

std::vector<std::string> SetgridHelp = {
  "NAME",
  "    setgrid, setgridtype, setgridarea - Set grid information",
  "",
  "SYNOPSIS",
  "    setgrid,grid  infile outfile",
  "    setgridtype,gridtype  infile outfile",
  "    setgridarea,gridarea  infile outfile",
  "",
  "DESCRIPTION",
  "    This module modifies the metadata of the horizontal grid. Depending on "
  "the ",
  "    chosen operator a new grid description is set, the coordinates are "
  "converted",
  "    or the grid cell area is added.",
  "",
  "OPERATORS",
  "    setgrid      Set grid",
  "                 Sets a new grid description. The input fields need to have "
  "the same grid size",
  "                 as the size of the target grid description.",
  "    setgridtype  Set grid type",
  "                 Sets the grid type of all input fields. The following grid "
  "types are available:",
  "                 curvilinear "
  "    Converts a regular grid to a curvilinear grid",
  "                 unstructured"
  "    Converts a regular or curvilinear grid to an unstructured grid",
  "                 dereference "
  "    Dereference a reference to a grid",
  "                 regular     "
  "    Linear interpolation of a reduced Gaussian grid to a regular Gaussian "
  "grid",
  "                 regularnn   "
  "    Nearest neighbor interpolation of a reduced Gaussian grid to a regular "
  "Gaussian grid",
  "                 lonlat      "
  "    Converts a regular lonlat grid stored as a curvilinear grid back to a "
  "lonlat grid",
  "    setgridarea  Set grid cell area",
  "                 Sets the grid cell area. The parameter gridarea is the "
  "path to a data file,",
  "                 the first field is used as grid cell area. The input "
  "fields need to have the same",
  "                 grid size as the grid cell area. The grid cell area is "
  "used to compute",
  "                 the weights of each grid cell if needed by an operator, "
  "e.g. for fldmean.",
  "",
  "PARAMETER",
  "    grid      STRING  Grid description file or name",
  "    gridtype  STRING  Grid type (curvilinear, unstructured, regular, lonlat "
  "or dereference)",
  "    gridarea  STRING  Data file, the first field is used as grid cell area",
};

std::vector<std::string> SetzaxisHelp = {
  "NAME",
  "    setzaxis, genlevelbounds - Set z-axis information",
  "",
  "SYNOPSIS",
  "    setzaxis,zaxis  infile outfile",
  "    genlevelbounds[,zbot[,ztop]]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module modifies the metadata of the vertical grid.",
  "",
  "OPERATORS",
  "    setzaxis        Set z-axis",
  "                    This operator sets the z-axis description of all "
  "variables with the same number of level as the new z-axis.",
  "    genlevelbounds  Generate level bounds",
  "                    Generates the layer bounds of the z-axis.",
  "",
  "PARAMETER",
  "    zaxis  STRING  Z-axis description file or name of the target z-axis",
  "    zbot   FLOAT   Specifying the bottom of the vertical column. Must have "
  "the same units as z-axis. ",
  "    ztop   FLOAT   Specifying the top of the vertical column. Must have the "
  "same units as z-axis. ",
};

std::vector<std::string> InvertHelp = {
  "NAME",
  "    invertlat - Invert latitudes",
  "",
  "SYNOPSIS",
  "    invertlat  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator inverts the latitudes of all fields on a rectilinear "
  "grid. ",
};

std::vector<std::string> InvertlevHelp = {
  "NAME",
  "    invertlev - Invert levels",
  "",
  "SYNOPSIS",
  "    invertlev  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator inverts the levels of all 3D variables.",
};

std::vector<std::string> ShiftxyHelp = {
  "NAME",
  "    shiftx, shifty - Shift field",
  "",
  "SYNOPSIS",
  "    <operator>,<nshift>,<cyclic>,<coord>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators to shift all fields in x or y direction.",
  "    All fields need to have the same horizontal rectilinear or curvilinear "
  "grid.",
  "",
  "OPERATORS",
  "    shiftx  Shift x",
  "            Shifts all fields in x direction.",
  "    shifty  Shift y",
  "            Shifts all fields in y direction.",
  "",
  "PARAMETER",
  "    nshift  INTEGER  Number of grid cells to shift (default: 1)",
  "    cyclic  STRING   If set, cells are filled up cyclic (default: missing "
  "value)",
  "    coord   STRING   If set, coordinates are also shifted",
};

std::vector<std::string> MaskregionHelp = {
  "NAME",
  "    maskregion - Mask regions",
  "",
  "SYNOPSIS",
  "    maskregion,regions  infile outfile",
  "",
  "DESCRIPTION",
  "    Masks different regions of fields with a regular lon/lat grid. The "
  "elements ",
  "    inside a region are untouched, the elements outside are set to missing "
  "value.",
  "    Considered are only those grid cells with the grid center inside the "
  "regions.",
  "    All input fields must have the same horizontal grid.",
  "    The user has to give ASCII formatted files with different regions.",
  "    A region is defined by a polygon. Each line of a polygon description "
  "file ",
  "    contains the longitude and latitude of one point.",
  "    Each polygon description file can contain one or more polygons "
  "separated",
  "    by a line with the character \\&.",
  "",
  "PARAMETER",
  "    regions  STRING Comma separated list of ASCII formatted files with "
  "different regions",
};

std::vector<std::string> MaskboxHelp = {
  "NAME",
  "    masklonlatbox, maskindexbox - Mask a box",
  "",
  "SYNOPSIS",
  "    masklonlatbox,lon1,lon2,lat1,lat2  infile outfile",
  "    maskindexbox,idx1,idx2,idy1,idy2  infile outfile",
  "",
  "DESCRIPTION",
  "    Masked a box of the rectangularly understood field. The elements inside "
  "the box are untouched, the ",
  "    elements outside are set to missing value. All input fields need to "
  "have the same horizontal grid.",
  "    Use sellonlatbox or selindexbox if only the data inside the box are "
  "needed.",
  "",
  "OPERATORS",
  "    masklonlatbox  Mask a longitude/latitude box",
  "                   Masked a regular longitude/latitude box. The user has to "
  "give the longitudes and latitudes of the ",
  "                   edges of the box. Considered are only those grid cells "
  "with the grid center inside the lon/lat box.",
  "    maskindexbox   Mask an index box",
  "                   Masked an index box. The user has to give the indexes of "
  "the edges of the box. ",
  "                   The index of the left edge can be greater then the one "
  "of the right edge.",
  "",
  "PARAMETER",
  "    lon1  FLOAT    Western longitude",
  "    lon2  FLOAT    Eastern longitude",
  "    lat1  FLOAT    Southern or northern latitude",
  "    lat2  FLOAT    Northern or southern latitude",
  "    idx1  INTEGER  Index of first longitude",
  "    idx2  INTEGER  Index of last longitude",
  "    idy1  INTEGER  Index of first latitude",
  "    idy2  INTEGER  Index of last latitude",
};

std::vector<std::string> SetboxHelp = {
  "NAME",
  "    setclonlatbox, setcindexbox - Set a box to constant",
  "",
  "SYNOPSIS",
  "    setclonlatbox,c,lon1,lon2,lat1,lat2  infile outfile",
  "    setcindexbox,c,idx1,idx2,idy1,idy2  infile outfile",
  "",
  "DESCRIPTION",
  "    Sets a box of the rectangularly understood field to a constant value. "
  "The elements outside ",
  "    the box are untouched, the elements inside are set to the given "
  "constant. All input fields ",
  "    need to have the same horizontal grid.",
  "",
  "OPERATORS",
  "    setclonlatbox  Set a longitude/latitude box to constant",
  "                   Sets the values of a longitude/latitude box to a "
  "constant value. The ",
  "                   user has to give the longitudes and latitudes of the "
  "edges of the box.",
  "    setcindexbox   Set an index box to constant",
  "                   Sets the values of an index box to a constant value. The "
  "user has to ",
  "                   give the indexes of the edges of the box. The index of "
  "the left edge ",
  "                   can be greater than the one of the right edge.",
  "",
  "PARAMETER",
  "    c     FLOAT    Constant",
  "    lon1  FLOAT    Western longitude",
  "    lon2  FLOAT    Eastern longitude",
  "    lat1  FLOAT    Southern or northern latitude",
  "    lat2  FLOAT    Northern or southern latitude",
  "    idx1  INTEGER  Index of first longitude",
  "    idx2  INTEGER  Index of last longitude",
  "    idy1  INTEGER  Index of first latitude",
  "    idy2  INTEGER  Index of last latitude",
};

std::vector<std::string> EnlargeHelp = {
  "NAME",
  "    enlarge - Enlarge fields",
  "",
  "SYNOPSIS",
  "    enlarge,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    Enlarge all fields of infile to a user given horizontal grid. Normally "
  "only the last ",
  "    field element is used for the enlargement. If however the input and "
  "output",
  "    grid are regular lon/lat grids, a zonal or meridional enlargement is "
  "possible.",
  "    Zonal enlargement takes place, if the xsize of the input field is 1 "
  "and ",
  "    the ysize of both grids are the same. For meridional enlargement the "
  "ysize",
  "    have to be 1 and the xsize of both grids should have the same size.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
};

std::vector<std::string> SetmissHelp = {
  "NAME",
  "    setmissval, setctomiss, setmisstoc, setrtomiss, setvrange, "
  "setmisstonn, ",
  "    setmisstodis - Set missing value",
  "",
  "SYNOPSIS",
  "    setmissval,newmiss  infile outfile",
  "    setctomiss,c  infile outfile",
  "    setmisstoc,c  infile outfile",
  "    setrtomiss,rmin,rmax  infile outfile",
  "    setvrange,rmin,rmax  infile outfile",
  "    setmisstonn  infile outfile",
  "    setmisstodis[,neighbors]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module sets part of a field to missing value or missing values",
  "    to a constant value. Which part of the field is set depends on the ",
  "    chosen operator.",
  "",
  "OPERATORS",
  "    setmissval    Set a new missing value",
  "                           / newmiss   if i(t,x) EQ miss",
  "                  o(t,x) = ",
  "                           \\ i(t,x)    if i(t,x) NE miss",
  "    setctomiss    Set constant to missing value",
  "                           / miss   if i(t,x) EQ c",
  "                  o(t,x) = ",
  "                           \\ i(t,x) if i(t,x) NE c",
  "    setmisstoc    Set missing value to constant",
  "                           / c      if i(t,x) EQ miss",
  "                  o(t,x) = ",
  "                           \\ i(t,x) if i(t,x) NE miss",
  "    setrtomiss    Set range to missing value",
  "                           / miss   if i(t,x) GE rmin AND i(t,x) LE rmax",
  "                  o(t,x) = ",
  "                           \\ i(t,x) if i(t,x) LT rmin OR  i(t,x) GT rmax",
  "    setvrange     Set valid range",
  "                           / miss   if i(t,x) LT rmin OR  i(t,x) GT rmax",
  "                  o(t,x) = ",
  "                           \\ i(t,x) if i(t,x) GE rmin AND i(t,x) LE rmax",
  "    setmisstonn   Set missing value to nearest neighbor",
  "                  Set all missing values to the nearest non missing value.",
  "                           / i(t,y) if i(t,x) EQ miss AND i(t,y) NE miss",
  "                  o(t,x) = ",
  "                           \\ i(t,x) if i(t,x) NE miss",
  "    setmisstodis  Set missing value to distance-weighted average",
  "                  Set all missing values to the distance-weighted average "
  "of the nearest non missing values.",
  "                  The default number of nearest neighbors is 4.",
  "",
  "PARAMETER",
  "    neighbors  INTEGER  Number of nearest neighbors",
  "    newmiss    FLOAT    New missing value",
  "    c          FLOAT    Constant",
  "    rmin       FLOAT    Lower bound",
  "    rmax       FLOAT    Upper bound",
};

std::vector<std::string> ExprHelp = {
  "NAME",
  "    expr, exprf, aexpr, aexprf - Evaluate expressions",
  "",
  "SYNOPSIS",
  "    expr,instr  infile outfile",
  "    exprf,filename  infile outfile",
  "    aexpr,instr  infile outfile",
  "    aexprf,filename  infile outfile",
  "",
  "DESCRIPTION",
  "    This module arithmetically processes every timestep of the input "
  "dataset.",
  "    Each individual assignment statement have to end with a semi-colon.",
  "    Unlike regular variables, temporary variables are never written to the "
  "output stream.",
  "    To define a temporary variable simply prefix the variable name with an "
  "underscore (e.g. _varname)",
  "    when the variable is declared.",
  "    ",
  "    The following operators are supported:",
  "    ",
  "     Operator   & Meaning             & Example   & Result ",
  "         =      & assignment          & x = y     & Assigns y to x",
  "         +      & addition            & x + y     & Sum of x and y",
  "         -      & subtraction         & x - y     & Difference of x and y   "
  " ",
  "         *      & multiplication      & x * y     & Product of x and y ",
  "         /      & division            & x / y     & Quotient of x and y",
  "         ^      & exponentiation      & x ^y      & Exponentiates x with y ",
  "         ==     & equal to            & x == y    &  1, if x equal to y; "
  "else 0",
  "         !=     & not equal to        & x != y    &  1, if x not equal to "
  "y; else 0",
  "         >      & greater than        & x > y     &  1, if x greater than "
  "y; else 0",
  "         <      & less than           & x < y     &  1, if x less than y; "
  "else 0",
  "         >=     & greater equal       & x >= y    &  1, if x greater equal "
  "y; else 0",
  "         <=     & less equal          & x <= y    &  1, if x less equal y; "
  "else 0",
  "         <=>    & less equal greater  & x <=> y   & -1, if x less y; 1, if "
  "x greater y; else 0 ",
  "         &&     & logical AND         & x && y    &  1, if x and y not "
  "equal 0; else 0",
  "         ||     & logical OR          & x || y    &  1, if x or y not equal "
  "0; else 0",
  "         !      & logical NOT         & !x        &  1, if x equal 0; else "
  "0",
  "         ?:     & ternary conditional & x ? y : z & y, if x not equal 0, "
  "else z ",
  "    ",
  "    The following functions are supported:",
  "    ",
  "    Math intrinsics:",
  "    ",
  "    abs(x)  "
  "    Absolute value of x",
  "    floor(x)"
  "    Round to largest integral value not greater than x",
  "    ceil(x) "
  "    Round to smallest integral value not less than x",
  "    float(x)"
  "    32-bit float value of x",
  "    int(x)  "
  "    Integer value of x",
  "    nint(x) "
  "    Nearest integer value of x",
  "    sqr(x)  "
  "    Square of x",
  "    sqrt(x) "
  "    Square Root of x",
  "    exp(x)  "
  "    Exponential of x",
  "    ln(x)   "
  "    Natural logarithm of x",
  "    log10(x)"
  "    Base 10 logarithm of x",
  "    sin(x)  "
  "    Sine of x, where x is specified in radians",
  "    cos(x)  "
  "    Cosine of x, where x is specified in radians",
  "    tan(x)  "
  "    Tangent of x, where x is specified in radians",
  "    asin(x) "
  "    Arc-sine of x, where x is specified in radians",
  "    acos(x) "
  "    Arc-cosine of x, where x is specified in radians",
  "    atan(x) "
  "    Arc-tangent of x, where x is specified in radians",
  "    rad(x)  "
  "    Convert x from degrees to radians",
  "    deg(x)  "
  "    Convert x from radians to degrees",
  "    ",
  "    Coordinates:",
  "    ",
  "    clon(x)    "
  "    Longitude coordinate of x (available only if x has geographical "
  "coordinates) ",
  "    clat(x)    "
  "    Latitude coordinate of x (available only if x has geographical "
  "coordinates) ",
  "    gridarea(x)"
  "    Grid cell area of x (available only if x has geographical coordinates) ",
  "    clev(x)    "
  "    Level coordinate of x (0, if x is a 2D surface variable)",
  "    ctimestep()"
  "    Timestep number (1 to N)",
  "    cdate()    "
  "    Verification date as YYYYMMDD",
  "    ctime()    "
  "    Verification time as HHMMSS.millisecond",
  "    cdeltat()  "
  "    Difference between current and last timestep in seconds",
  "    cday()     "
  "    Day as DD",
  "    cmonth()   "
  "    Month as MM",
  "    cyear()    "
  "    Year as YYYY",
  "    csecond()  "
  "    Second as SS.millisecond",
  "    cminute()  "
  "    Minute as MM",
  "    chour()    "
  "    Hour as HH",
  "    ",
  "    Constants:",
  "    ",
  "    ngp(x)    "
  "    Number of horizontal grid points",
  "    nlev(x)   "
  "    Number of vertical levels",
  "    size(x)   "
  "    Total number of elements (ngp(x)*nlev(x))",
  "    missval(x)"
  "    Returns the missing value of variable x",
  "    ",
  "    Statistical values over a field:",
  "    ",
  "    fldmin(x), fldmax(x), fldsum(x), fldmean(x), fldavg(x), fldstd(x), "
  "fldstd1(x), fldvar(x), fldvar1(x)",
  "    ",
  "    Vertical statistical values:",
  "    ",
  "    vertmin(x), vertmax(x), vertsum(x), vertmean(x), vertavg(x), "
  "vertstd(x), vertstd1(x), vertvar(x), vertvar1(x)",
  "    ",
  "    Miscellaneous:",
  "    ",
  "    sellevel(x,k) "
  "    Select level k of variable x",
  "    sellevidx(x,k)"
  "    Select level index k of variable x",
  "    remove(x)     "
  "    Remove variable x from output stream",
  "    ",
  "",
  "OPERATORS",
  "    expr    Evaluate expressions",
  "            The processing instructions are read from the parameter.",
  "    exprf   Evaluate expressions script",
  "            Contrary to expr the processing instructions are read from a "
  "file.",
  "    aexpr   Evaluate expressions and append results",
  "            Same as expr, but keep input variables and append results",
  "    aexprf  Evaluate expression script and append results",
  "            Same as exprf, but keep input variables and append results",
  "",
  "PARAMETER",
  "    instr     STRING  Processing instructions (need to be 'quoted' in most "
  "cases)",
  "    filename  STRING  File with processing instructions",
  "",
  "NOTE",
  "    The expr commands sellevel(x,k) and sellevidx(x,k) are only available "
  "with exprf/aexprf.",
  "    If the input stream contains duplicate entries of the same variable "
  "name then the last one is used.",
};

std::vector<std::string> MathHelp = {
  "NAME",
  "    abs, int, nint, pow, sqr, sqrt, exp, ln, log10, sin, cos, tan, asin, "
  "acos, ",
  "    atan, reci, not - Mathematical functions",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains some standard mathematical functions.",
  "    All trigonometric functions calculate with radians.",
  "",
  "OPERATORS",
  "    abs    Absolute value",
  "           o(t,x) = abs(i(t,x))",
  "    int    Integer value",
  "           o(t,x) = int(i(t,x))",
  "    nint   Nearest integer value",
  "           o(t,x) = nint(i(t,x))",
  "    pow    Power",
  "           o(t,x) = i(t,x)^y",
  "    sqr    Square",
  "           o(t,x) = i(t,x)^2",
  "    sqrt   Square root",
  "           o(t,x) = sqrt(i(t,x))",
  "    exp    Exponential",
  "           o(t,x) = e^i(t,x)",
  "    ln     Natural logarithm",
  "           o(t,x) = ln(i(t,x))",
  "    log10  Base 10 logarithm",
  "           o(t,x) = log10(i(t,x))",
  "    sin    Sine",
  "           o(t,x) = sin(i(t,x))",
  "    cos    Cosine",
  "           o(t,x) = cos(i(t,x))",
  "    tan    Tangent",
  "           o(t,x) = tan(i(t,x))",
  "    asin   Arc sine",
  "           o(t,x) = asin(i(t,x))",
  "    acos   Arc cosine",
  "           o(t,x) = acos(i(t,x))",
  "    atan   Arc tangent",
  "           o(t,x) = atan(i(t,x))",
  "    reci   Reciprocal value",
  "           o(t,x) = 1 / i(t,x)",
  "    not    Logical NOT",
  "           o(t,x) = 1, if x equal 0; else 0",
};

std::vector<std::string> ArithcHelp = {
  "NAME",
  "    addc, subc, mulc, divc - Arithmetic with a constant",
  "",
  "SYNOPSIS",
  "    <operator>,c  infile outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic with all field elements of a "
  "dataset and ",
  "    a constant. The fields in outfile inherit the meta data from infile.",
  "",
  "OPERATORS",
  "    addc  Add a constant",
  "          o(t,x) = i(t,x) + c",
  "    subc  Subtract a constant",
  "          o(t,x) = i(t,x) - c",
  "    mulc  Multiply with a constant",
  "          o(t,x) = i(t,x) * c",
  "    divc  Divide by a constant",
  "          o(t,x) = i(t,x) / c",
  "",
  "PARAMETER",
  "    c  FLOAT  Constant",
};

std::vector<std::string> ArithHelp = {
  "NAME",
  "    add, sub, mul, div, min, max, atan2 - Arithmetic on two datasets",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic of two datasets.",
  "    The number of fields in infile1 should be the same as in infile2.",
  "    The fields in outfile inherit the meta data from infile1.",
  "    One of the input files can contain only one timestep or one variable.",
  "",
  "OPERATORS",
  "    add    Add two fields",
  "           o(t,x) = i_1(t,x) + i_2(t,x)",
  "    sub    Subtract two fields",
  "           o(t,x) = i_1(t,x) - i_2(t,x)",
  "    mul    Multiply two fields",
  "           o(t,x) = i_1(t,x) * i_2(t,x)",
  "    div    Divide two fields",
  "           o(t,x) = i_1(t,x) / i_2(t,x)",
  "    min    Minimum of two fields",
  "           o(t,x) = min(i_1(t,x), i_2(t,x))",
  "    max    Maximum of two fields",
  "           o(t,x) = max(i_1(t,x), i_2(t,x))",
  "    atan2  Arc tangent of two fields",
  "           The atan2 operator calculates the arc tangent of two fields. The "
  "result is",
  "           in radians, which is between -PI and PI (inclusive).",
  "           ",
  "           o(t,x) = atan2(i_1(t,x), i_2(t,x))",
};

std::vector<std::string> MonarithHelp = {
  "NAME",
  "    monadd, monsub, monmul, mondiv - Monthly arithmetic",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic of a time series and one",
  "    timestep with the same month and year. For each field in infile1",
  "    the corresponding field of the timestep in infile2 with the",
  "    same month and year is used. The header information in infile1",
  "    have to be the same as in infile2. Usually infile2 is generated",
  "    by an operator of the module MONSTAT.",
  "",
  "OPERATORS",
  "    monadd  Add monthly time series",
  "            Adds a time series and a monthly time series.",
  "    monsub  Subtract monthly time series",
  "            Subtracts a time series and a monthly time series.",
  "    monmul  Multiply monthly time series",
  "            Multiplies a time series and a monthly time series.",
  "    mondiv  Divide monthly time series",
  "            Divides a time series and a monthly time series.",
};

std::vector<std::string> YhourarithHelp = {
  "NAME",
  "    yhouradd, yhoursub, yhourmul, yhourdiv - Multi-year hourly arithmetic",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic of a time series and one",
  "    timestep with the same hour and day of year. For each field in infile1",
  "    the corresponding field of the timestep in infile2 with the",
  "    same hour and day of year is used. The header information in infile1",
  "    have to be the same as in infile2. Usually infile2 is generated",
  "    by an operator of the module YHOURSTAT.",
  "",
  "OPERATORS",
  "    yhouradd  Add multi-year hourly time series",
  "              Adds a time series and a multi-year hourly time series.",
  "    yhoursub  Subtract multi-year hourly time series",
  "              Subtracts a time series and a multi-year hourly time series.",
  "    yhourmul  Multiply multi-year hourly time series",
  "              Multiplies a time series and a multi-year hourly time series.",
  "    yhourdiv  Divide multi-year hourly time series",
  "              Divides a time series and a multi-year hourly time series.",
};

std::vector<std::string> YdayarithHelp = {
  "NAME",
  "    ydayadd, ydaysub, ydaymul, ydaydiv - Multi-year daily arithmetic",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic of a time series and one",
  "    timestep with the same day of year. For each field in infile1",
  "    the corresponding field of the timestep in infile2 with the",
  "    same day of year is used. The header information in infile1",
  "    have to be the same as in infile2. Usually infile2 is generated",
  "    by an operator of the module YDAYSTAT.",
  "",
  "OPERATORS",
  "    ydayadd  Add multi-year daily time series",
  "             Adds a time series and a multi-year daily time series.",
  "    ydaysub  Subtract multi-year daily time series",
  "             Subtracts a time series and a multi-year daily time series.",
  "    ydaymul  Multiply multi-year daily time series",
  "             Multiplies a time series and a multi-year daily time series.",
  "    ydaydiv  Divide multi-year daily time series",
  "             Divides a time series and a multi-year daily time series.",
};

std::vector<std::string> YmonarithHelp = {
  "NAME",
  "    ymonadd, ymonsub, ymonmul, ymondiv - Multi-year monthly arithmetic",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic of a time series and one "
  "timestep",
  "    with the same month of year. For each field in infile1 the "
  "corresponding",
  "    field of the timestep in infile2 with the same month of year is used.",
  "    The header information in infile1 have to be the same as in infile2.",
  "    Usually infile2 is generated by an operator of the module YMONSTAT.",
  "",
  "OPERATORS",
  "    ymonadd  Add multi-year monthly time series",
  "             Adds a time series and a multi-year monthly time series.",
  "    ymonsub  Subtract multi-year monthly time series",
  "             Subtracts a time series and a multi-year monthly time series.",
  "    ymonmul  Multiply multi-year monthly time series",
  "             Multiplies a time series and a multi-year monthly time series.",
  "    ymondiv  Divide multi-year monthly time series",
  "             Divides a time series and a multi-year monthly time series.",
};

std::vector<std::string> YseasarithHelp = {
  "NAME",
  "    yseasadd, yseassub, yseasmul, yseasdiv - Multi-year seasonal arithmetic",
  "",
  "SYNOPSIS",
  "    <operator>  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This module performs simple arithmetic of a time series and one "
  "timestep",
  "    with the same season. For each field in infile1 the corresponding",
  "    field of the timestep in infile2 with the same season is used.",
  "    The header information in infile1 have to be the same as in infile2.",
  "    Usually infile2 is generated by an operator of the module YSEASSTAT.",
  "",
  "OPERATORS",
  "    yseasadd  Add multi-year seasonal time series",
  "              Adds a time series and a multi-year seasonal time series.",
  "    yseassub  Subtract multi-year seasonal time series",
  "              Subtracts a time series and a multi-year seasonal time "
  "series.",
  "    yseasmul  Multiply multi-year seasonal time series",
  "              Multiplies a time series and a multi-year seasonal time "
  "series.",
  "    yseasdiv  Divide multi-year seasonal time series",
  "              Divides a time series and a multi-year seasonal time series.",
};

std::vector<std::string> ArithdaysHelp = {
  "NAME",
  "    muldpm, divdpm, muldpy, divdpy - Arithmetic with days",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module multiplies or divides each timestep of a dataset with the "
  "corresponding",
  "    days per month or days per year. The result of these functions depends "
  "on the used",
  "    calendar of the input data.",
  "",
  "OPERATORS",
  "    muldpm  Multiply with days per month",
  "            o(t,x) = i(t,x) * days_per_month",
  "    divdpm  Divide by days per month",
  "            o(t,x) = i(t,x) / days_per_month",
  "    muldpy  Multiply with days per year",
  "            o(t,x) = i(t,x) * days_per_year",
  "    divdpy  Divide by days per year",
  "            o(t,x) = i(t,x) / days_per_year",
};

std::vector<std::string> TimcumsumHelp = {
  "NAME",
  "    timcumsum - Cumulative sum over all timesteps",
  "",
  "SYNOPSIS",
  "    timcumsum  infile outfile",
  "",
  "DESCRIPTION",
  "    The timcumsum operator calculates the cumulative sum over all "
  "timesteps.",
  "    Missing values are treated as numeric zero when summing.",
  "    ",
  "    o(t,x) = sum{i(t',x), 0<t'<=t}",
};

std::vector<std::string> ConsecstatHelp = {
  "NAME",
  "    consecsum, consects - Consecute timestep periods",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes periods over all timesteps in infile where a",
  "    certain property is valid. The property can be chosen by creating a "
  "mask from",
  "    the original data, which is the expected input format for operators of "
  "this",
  "    module. Depending on the operator full information about each period or",
  "    just its length and ending date are computed.",
  "",
  "OPERATORS",
  "    consecsum  Consecutive Sum",
  "               This operator computes periods of consecutive timesteps "
  "similar to a",
  "               runsum, but periods are finished, when the mask value is 0. "
  "That way",
  "               multiple periods can be found. Timesteps from the input are "
  "preserved. Missing",
  "               values are handled like 0, i.e. finish periods of "
  "consecutive timesteps.",
  "    consects   Consecutive Timesteps",
  "               In contrast to the operator above consects only computes the "
  "length of each",
  "               period together with its last timestep. To be able to "
  "perform statistical",
  "               analysis like min, max or mean, everything else is set to "
  "missing value.",
};

std::vector<std::string> EnsstatHelp = {
  "NAME",
  "    ensmin, ensmax, ensrange, enssum, ensmean, ensavg, ensstd, ensstd1, "
  "ensvar, ",
  "    ensvar1, enspctl - Statistical values over an ensemble",
  "",
  "SYNOPSIS",
  "    <operator>  infiles outfile",
  "    enspctl,p  infiles outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over an ensemble of input "
  "files.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance,",
  "    standard deviation or a certain percentile over all input files is "
  "written",
  "    to outfile.",
  "    All input files need to have the same structure with the same "
  "variables.",
  "    The date information of a timestep in outfile is the date of the first "
  "input file.",
  "",
  "OPERATORS",
  "    ensmin    Ensemble minimum",
  "              o(t,x) = min{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensmax    Ensemble maximum",
  "              o(t,x) = max{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensrange  Ensemble range",
  "              o(t,x) = range{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    enssum    Ensemble sum",
  "              o(t,x) = sum{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensmean   Ensemble mean",
  "              o(t,x) = mean{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensavg    Ensemble average",
  "              o(t,x) = avg{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensstd    Ensemble standard deviation",
  "              Normalize by n.",
  "              ",
  "              o(t,x) = std{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensstd1   Ensemble standard deviation (n-1)",
  "              Normalize by (n-1).",
  "              ",
  "              o(t,x) = std1{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensvar    Ensemble variance",
  "              Normalize by n.",
  "              ",
  "              o(t,x) = var{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    ensvar1   Ensemble variance (n-1)",
  "              Normalize by (n-1).",
  "              ",
  "              o(t,x) = var1{i1(t,x), i2(t,x), ..., in(t,x)}",
  "    enspctl   Ensemble percentiles",
  "              o(t,x) = pth percentile {i1(t,x), i2(t,x), ..., in(t,x)}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "NOTE",
  "    This operator needs to open all input files simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> Ensstat2Help = {
  "NAME",
  "    ensrkhistspace, ensrkhisttime, ensroc - Statistical values over an "
  "ensemble",
  "",
  "SYNOPSIS",
  "    <operator>  obsfile ensfiles outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over the ensemble of ensfiles "
  "using",
  "    obsfile as a reference. Depending on the operator a ranked Histogram "
  "or ",
  "    a roc-curve over all Ensembles ensfiles",
  "    with reference to obsfile is written to outfile. ",
  "    The date and grid information of a timestep in outfile is the date of "
  "the ",
  "    first input file. Thus all input files are required to have the same "
  "structure in ",
  "    terms of the gridsize, variable definitions and number of timesteps. ",
  "    ",
  "    All Operators in this module use obsfile as the reference (for "
  "instance ",
  "    an observation) whereas ensfiles are understood as an ensemble "
  "consisting ",
  "    of n (where n is the number of ensfiles) members. ",
  "    ",
  "    The operators ensrkhistspace and ensrkhisttime compute Ranked "
  "Histograms. ",
  "    Therefor the vertical axis is utilized as the Histogram axis, which "
  "prohibits",
  "    the use of files containing more than one level. The histogram axis "
  "has ",
  "    nensfiles+1 bins with level 0 containing for each grid point the number "
  "of ",
  "    observations being smaller as all ensembles and level nensfiles+1 "
  "indicating",
  "    the number of observations being larger than all ensembles. ",
  "    ",
  "    ensrkhistspace computes a ranked histogram at each timestep reducing "
  "each ",
  "    horizontal grid to a 1x1 grid and keeping the time axis as in obsfile. ",
  "    Contrary ensrkhistspace  computes a histogram at each grid point "
  "keeping the ",
  "    horizontal grid for each variable and reducing the time-axis. The time "
  "information",
  "    is that from the last timestep in obsfile. ",
  "",
  "OPERATORS",
  "    ensrkhistspace  Ranked Histogram averaged over time",
  "    ensrkhisttime   Ranked Histogram averaged over space",
  "    ensroc          Ensemble Receiver Operating characteristics",
};

std::vector<std::string> EnsvalHelp = {
  "NAME",
  "    enscrps, ensbrs - Ensemble validation tools",
  "",
  "SYNOPSIS",
  "    enscrps  rfile infiles outfilebase",
  "    ensbrs,x  rfile infiles outfilebase",
  "",
  "DESCRIPTION",
  "    This module computes ensemble validation scores and their decomposition "
  "such as ",
  "    the Brier and cumulative ranked probability score (CRPS). ",
  "    The first file is used as a reference it can be a climatology, "
  "observation or ",
  "    reanalysis against which the skill of the ensembles given in infiles is "
  "measured. ",
  "    Depending on the operator a number of output files is generated each "
  "containing ",
  "    the skill score and its decomposition corresponding to the operator. ",
  "    The output is averaged  over horizontal fields using appropriate "
  "weights ",
  "    for each level and timestep in rfile. ",
  "    ",
  "    All input files need to have the same structure with the same "
  "variables.",
  "    The date information of a timestep in outfile is the date of the first "
  "input file.",
  "    The output files are named as ",
  "    <outfilebase>.<type>.<filesuffix> where <type> depends on the ",
  "    operator and <filesuffix> is determined from the output file type. ",
  "    There are three output files for operator enscrps and four output "
  "files ",
  "    for operator ensbrs.",
  "    ",
  "    The CRPS and its decomposition into Reliability and the potential ",
  "    CRPS are calculated by an appropriate averaging over the field ",
  "    members (note, that the CRPS does *not* average linearly). ",
  "    In the three output files ",
  "    <type> has the following meaning:",
  "    crps for the CRPS, reli for the reliability ",
  "    and crpspot for the potential crps. The relation ",
  "    CRPS = CRPS_{pot} + RELI",
  "    holds. 	  ",
  "    ",
  "    The Brier score of the Ensemble given by infiles with respect to the ",
  "    reference given in rfile and the threshold x is calculated. ",
  "    In the four output files <type> has the following meaning: ",
  "    brs for the Brier score wrt threshold  x; ",
  "    brsreli for the Brier score reliability wrt threshold x;",
  "    brsreso for the Brier score resolution wrt threshold x;",
  "    brsunct for the Brier score uncertainty wrt threshold x.",
  "    In analogy to the CRPS the following relation holds:",
  "    BRS(x) = RELI(x)-RESO(x)+ UNCT(x).",
  "    ",
  "    The implementation of the decomposition of the CRPS and Brier Score "
  "follows ",
  "      Hans Hersbach (2000): Decomposition of the Continuous Ranked "
  "Probability ",
  "      Score for Ensemble Prediction Systems, in: Weather and Forecasting "
  "(15) ",
  "      pp. 559-570. ",
  "    ",
  "    The CRPS code decomposition has been verified against the CRAN - "
  "ensemble ",
  "    validation package from R. Differences occur when grid-cell area is "
  "not ",
  "    uniform as the implementation in R does not account for that. ",
  "    ",
  "",
  "OPERATORS",
  "    enscrps  Ensemble CRPS and decomposition",
  "    ensbrs   Ensemble Brier score",
  "             Ensemble Brier Score and Decomposition",
};

std::vector<std::string> FldstatHelp = {
  "NAME",
  "    fldmin, fldmax, fldrange, fldsum, fldmean, fldavg, fldstd, fldstd1, "
  "fldvar, ",
  "    fldvar1, fldpctl - Statistical values over a field",
  "",
  "SYNOPSIS",
  "    <operator>,weights  infile outfile",
  "    fldpctl,p  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values of the input fields. According "
  "to the chosen ",
  "    operator the field minimum, maximum, range, sum, average, variance, "
  "standard deviation or ",
  "    a certain percentile is written to outfile.",
  "",
  "OPERATORS",
  "    fldmin    Field minimum",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = min{i(t,x'), x_1<x'<=x_n}",
  "    fldmax    Field maximum",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = max{i(t,x'), x_1<x'<=x_n}",
  "    fldrange  Field range",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = range{i(t,x'), x_1<x'<=x_n}",
  "    fldsum    Field sum",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = sum{i(t,x'), x_1<x'<=x_n}",
  "    fldmean   Field mean",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = mean{i(t,x'), x_1<x'<=x_n}",
  "              weighted by area weights obtained by the input field.",
  "    fldavg    Field average",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = avg{i(t,x'), x_1<x'<=x_n}",
  "              weighted by area weights obtained by the input field.",
  "    fldstd    Field standard deviation",
  "              Normalize by n. For every gridpoint x_1, ..., x_n of the same "
  "field it is:",
  "              ",
  "              o(t,1) = std{i(t,x'), x_1<x'<=x_n}",
  "              weighted by area weights obtained by the input field.",
  "    fldstd1   Field standard deviation (n-1)",
  "              Normalize by (n-1). For every gridpoint x_1, ..., x_n of the "
  "same field it is:",
  "              ",
  "              o(t,1) = std1{i(t,x'), x_1<x'<=x_n}",
  "              weighted by area weights obtained by the input field.",
  "    fldvar    Field variance",
  "              Normalize by n. For every gridpoint x_1, ..., x_n of the same "
  "field it is:",
  "              ",
  "              o(t,1) = var{i(t,x'), x_1<x'<=x_n}",
  "              weighted by area weights obtained by the input field.",
  "    fldvar1   Field variance (n-1)",
  "              Normalize by (n-1). For every gridpoint x_1, ..., x_n of the "
  "same field it is:",
  "              ",
  "              o(t,1) = var1{i(t,x'), x_1<x'<=x_n}",
  "              weighted by area weights obtained by the input field.",
  "    fldpctl   Field percentiles",
  "              For every gridpoint x_1, ..., x_n of the same field it is:",
  "              ",
  "              o(t,1) = pth percentile {i(t,x'), x_1<x'<=x_n}",
  "",
  "PARAMETER",
  "    weights  BOOL   weights=FALSE disables weighting by grid cell area "
  "[default: weights=TRUE]",
  "    p        FLOAT  Percentile number in {0, ..., 100}",
};

std::vector<std::string> ZonstatHelp = {
  "NAME",
  "    zonmin, zonmax, zonrange, zonsum, zonmean, zonavg, zonstd, zonstd1, "
  "zonvar, ",
  "    zonvar1, zonpctl - Zonal statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "    zonpctl,p  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes zonal statistical values of the input fields.",
  "    According to the chosen operator the zonal minimum, maximum, range, "
  "sum, average,",
  "    variance, standard deviation or a certain percentile is written to "
  "outfile.",
  "    This operator requires all variables on the same regular lon/lat grid.",
  "",
  "OPERATORS",
  "    zonmin    Zonal minimum",
  "              For every latitude the minimum over all longitudes is "
  "computed.",
  "    zonmax    Zonal maximum",
  "              For every latitude the maximum over all longitudes is "
  "computed.",
  "    zonrange  Zonal range",
  "              For every latitude the range over all longitudes is computed.",
  "    zonsum    Zonal sum",
  "              For every latitude the sum over all longitudes is computed.",
  "    zonmean   Zonal mean",
  "              For every latitude the mean over all longitudes is computed.",
  "    zonavg    Zonal average",
  "              For every latitude the average over all longitudes is "
  "computed.",
  "    zonstd    Zonal standard deviation",
  "              For every latitude the standard deviation over all longitudes "
  "is computed. Normalize by n.",
  "    zonstd1   Zonal standard deviation (n-1)",
  "              For every latitude the standard deviation over all longitudes "
  "is computed. Normalize by (n-1). ",
  "    zonvar    Zonal variance",
  "              For every latitude the variance over all longitudes is "
  "computed. Normalize by n.",
  "    zonvar1   Zonal variance (n-1)",
  "              For every latitude the variance over all longitudes is "
  "computed. Normalize by (n-1).",
  "    zonpctl   Zonal percentiles",
  "              For every latitude the pth percentile over all longitudes is "
  "computed.",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
};

std::vector<std::string> MerstatHelp = {
  "NAME",
  "    mermin, mermax, merrange, mersum, mermean, meravg, merstd, merstd1, "
  "mervar, ",
  "    mervar1, merpctl - Meridional statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "    merpctl,p  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes meridional statistical values of the input fields.",
  "    According to the chosen operator the meridional minimum, maximum, "
  "range, sum, average,",
  "    variance, standard deviation or a certain percentile is written to "
  "outfile.",
  "    This operator requires all variables on the same regular lon/lat grid.",
  "",
  "OPERATORS",
  "    mermin    Meridional minimum",
  "              For every longitude the minimum over all latitudes is "
  "computed.",
  "    mermax    Meridional maximum",
  "              For every longitude the maximum over all latitudes is "
  "computed.",
  "    merrange  Meridional range",
  "              For every longitude the range over all latitudes is computed.",
  "    mersum    Meridional sum",
  "              For every longitude the sum over all latitudes is computed.",
  "    mermean   Meridional mean",
  "              For every longitude the area weighted mean over all latitudes "
  "is computed.",
  "    meravg    Meridional average",
  "              For every longitude the area weighted average over all "
  "latitudes is computed.",
  "    merstd    Meridional standard deviation",
  "              For every longitude the standard deviation over all latitudes "
  "is computed. Normalize by n.",
  "    merstd1   Meridional standard deviation (n-1)",
  "              For every longitude the standard deviation over all latitudes "
  "is computed. Normalize by (n-1).",
  "    mervar    Meridional variance",
  "              For every longitude the variance over all latitudes is "
  "computed. Normalize by n.",
  "    mervar1   Meridional variance (n-1)",
  "              For every longitude the variance over all latitudes is "
  "computed. Normalize by (n-1).",
  "    merpctl   Meridional percentiles",
  "              For every longitude the pth percentile over all latitudes is "
  "computed.",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
};

std::vector<std::string> GridboxstatHelp = {
  "NAME",
  "    gridboxmin, gridboxmax, gridboxrange, gridboxsum, gridboxmean, "
  "gridboxavg, ",
  "    gridboxstd, gridboxstd1, gridboxvar, gridboxvar1 - ",
  "    Statistical values over grid boxes",
  "",
  "SYNOPSIS",
  "    <operator>,nx,ny  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over surrounding grid boxes.",
  "    According to the chosen operator the minimum, maximum, range, sum, "
  "average, ",
  "    variance, or standard deviation of the neighboring grid boxes is "
  "written to outfile.",
  "    All gridbox operators only works on quadrilateral curvilinear grids.",
  "",
  "OPERATORS",
  "    gridboxmin    Gridbox minimum",
  "                  Minimum value of the selected grid boxes.",
  "    gridboxmax    Gridbox maximum",
  "                  Maximum value of the selected grid boxes.",
  "    gridboxrange  Gridbox range",
  "                  Range (max-min value) of the selected grid boxes.",
  "    gridboxsum    Gridbox sum",
  "                  Sum of the selected grid boxes.",
  "    gridboxmean   Gridbox mean",
  "                  Mean of the selected grid boxes.",
  "    gridboxavg    Gridbox average",
  "                  Average of the selected grid boxes.",
  "    gridboxstd    Gridbox standard deviation",
  "                  Standard deviation of the selected grid boxes. Normalize "
  "by n.",
  "    gridboxstd1   Gridbox standard deviation (n-1)",
  "                  Standard deviation of the selected grid boxes. Normalize "
  "by (n-1).",
  "    gridboxvar    Gridbox variance",
  "                  Variance of the selected grid boxes. Normalize by n.",
  "    gridboxvar1   Gridbox variance (n-1)",
  "                  Variance of the selected grid boxes. Normalize by (n-1).",
  "",
  "PARAMETER",
  "    nx  INTEGER  Number of grid boxes in x direction",
  "    ny  INTEGER  Number of grid boxes in y direction",
};

std::vector<std::string> VertstatHelp = {
  "NAME",
  "    vertmin, vertmax, vertrange, vertsum, vertmean, vertavg, vertstd, "
  "vertstd1, ",
  "    vertvar, vertvar1 - Vertical statistical values",
  "",
  "SYNOPSIS",
  "    <operator>,weights  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over all levels of the input "
  "variables.",
  "    According to chosen operator the vertical minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation is written to outfile.",
  "",
  "OPERATORS",
  "    vertmin    Vertical minimum",
  "               For every gridpoint the minimum over all levels is computed.",
  "    vertmax    Vertical maximum",
  "               For every gridpoint the maximum over all levels is computed.",
  "    vertrange  Vertical range",
  "               For every gridpoint the range over all levels is computed.",
  "    vertsum    Vertical sum",
  "               For every gridpoint the sum over all levels is computed.",
  "    vertmean   Vertical mean",
  "               For every gridpoint the layer weighted mean over all levels "
  "is computed.",
  "    vertavg    Vertical average",
  "               For every gridpoint the layer weighted average over all "
  "levels is computed.",
  "    vertstd    Vertical standard deviation",
  "               For every gridpoint the standard deviation over all levels "
  "is computed. Normalize by n.",
  "    vertstd1   Vertical standard deviation (n-1)",
  "               For every gridpoint the standard deviation over all levels "
  "is computed. Normalize by (n-1).",
  "    vertvar    Vertical variance",
  "               For every gridpoint the variance over all levels is "
  "computed. Normalize by n.",
  "    vertvar1   Vertical variance (n-1)",
  "               For every gridpoint the variance over all levels is "
  "computed. Normalize by (n-1).",
  "",
  "PARAMETER",
  "    weights  BOOL   weights=FALSE disables weighting by layer thickness "
  "[default: weights=TRUE]",
};

std::vector<std::string> TimselstatHelp = {
  "NAME",
  "    timselmin, timselmax, timselrange, timselsum, timselmean, timselavg, ",
  "    timselstd, timselstd1, timselvar, timselvar1 - Time range statistical "
  "values",
  "",
  "SYNOPSIS",
  "    <operator>,nsets[,noffset[,nskip]]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values for a selected number of "
  "timesteps. According to ",
  "    the chosen operator the minimum, maximum, range, sum, average, variance "
  "or standard deviation of ",
  "    the selected timesteps is written to outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    timselmin    Time selection minimum",
  "                 For every adjacent sequence t1, ...., tn of timesteps of "
  "the same selected time range it is:",
  "                 ",
  "                 o(t,x) = min{i(t',x), t1 < t' <= tn}",
  "    timselmax    Time selection maximum",
  "                 For every adjacent sequence t1, ...., tn of timesteps of "
  "the same selected time range it is:",
  "                 ",
  "                 o(t,x) = max{i(t',x), t1 < t' <= tn}",
  "    timselrange  Time selection range",
  "                 For every adjacent sequence t1, ...., tn of timesteps of "
  "the same selected time range it is:",
  "                 ",
  "                 o(t,x) = range{i(t',x), t1 < t' <= tn}",
  "    timselsum    Time selection sum",
  "                 For every adjacent sequence t1, ...., tn of timesteps of "
  "the same selected time range it is:",
  "                 ",
  "                 o(t,x) = sum{i(t',x), t1 < t' <= tn}",
  "    timselmean   Time selection mean",
  "                 For every adjacent sequence t1, ...., tn of timesteps of "
  "the same selected time range it is:",
  "                 ",
  "                 o(t,x) = mean{i(t',x), t1 < t' <= tn}",
  "    timselavg    Time selection average",
  "                 For every adjacent sequence t1, ...., tn of timesteps of "
  "the same selected time range it is:",
  "                 ",
  "                 o(t,x) = avg{i(t',x), t1 < t' <= tn}",
  "    timselstd    Time selection standard deviation",
  "                 Normalize by n. For every adjacent sequence t1, ...., tn "
  "of timesteps of the same selected time range it is:",
  "                 ",
  "                 o(t,x) = std{i(t',x), t1 < t' <= tn}",
  "    timselstd1   Time selection standard deviation (n-1)",
  "                 Normalize by (n-1). For every adjacent sequence t1, ...., "
  "tn of timesteps of the same selected time range it is:",
  "                 ",
  "                 o(t,x) = std1{i(t',x), t1 < t' <= tn}",
  "    timselvar    Time selection variance",
  "                 Normalize by n. For every adjacent sequence t1, ...., tn "
  "of timesteps of the same selected time range it is:",
  "                 ",
  "                 o(t,x) = var{i(t',x), t1 < t' <= tn}",
  "    timselvar1   Time selection variance (n-1)",
  "                 Normalize by (n-1). For every adjacent sequence t1, ...., "
  "tn of timesteps of the same selected time range it is:",
  "                 ",
  "                 o(t,x) = var1{i(t',x), t1 < t' <= tn}",
  "",
  "PARAMETER",
  "    nsets    INTEGER  Number of input timesteps for each output timestep ",
  "    noffset  INTEGER  Number of input timesteps skipped before the first "
  "timestep range (optional)",
  "    nskip    INTEGER  Number of input timesteps skipped between timestep "
  "ranges (optional)",
};

std::vector<std::string> TimselpctlHelp = {
  "NAME",
  "    timselpctl - Time range percentile values",
  "",
  "SYNOPSIS",
  "    timselpctl,p,nsets[,noffset[,nskip]]  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentile values over a selected number of "
  "timesteps in infile1.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and infile3,",
  "    respectively. The default number of histogram bins is 101. The default "
  "can be overridden by setting the",
  "    environment variable CDO_PCTL_NBINS to a different value. The files "
  "infile2 and infile3 ",
  "    should be the result of corresponding timselmin and timselmax "
  "operations, respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    For every adjacent sequence t1, ...., tn of timesteps of the same "
  "selected time range it is:",
  "    ",
  "    o(t,x) = pth percentile {i(t',x), t1 < t' <= tn}",
  "",
  "PARAMETER",
  "    p        FLOAT    Percentile number in {0, ..., 100}",
  "    nsets    INTEGER  Number of input timesteps for each output timestep ",
  "    noffset  INTEGER  Number of input timesteps skipped before the first "
  "timestep range (optional)",
  "    nskip    INTEGER  Number of input timesteps skipped between timestep "
  "ranges (optional)",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> RunstatHelp = {
  "NAME",
  "    runmin, runmax, runrange, runsum, runmean, runavg, runstd, runstd1, "
  "runvar, ",
  "    runvar1 - Running statistical values",
  "",
  "SYNOPSIS",
  "    <operator>,nts  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes running statistical values over a selected number "
  "of timesteps. Depending on ",
  "    the chosen operator the minimum, maximum, range, sum, average, variance "
  "or standard deviation of a selected ",
  "    number of consecutive timesteps read from infile is written to "
  "outfile. ",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    runmin    Running minimum",
  "              o(t+(nts-1)/2,x) = min{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runmax    Running maximum",
  "              o(t+(nts-1)/2,x) = max{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runrange  Running range",
  "              o(t+(nts-1)/2,x) = range{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runsum    Running sum",
  "              o(t+(nts-1)/2,x) = sum{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runmean   Running mean",
  "              o(t+(nts-1)/2,x) = mean{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runavg    Running average",
  "              o(t+(nts-1)/2,x) = avg{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runstd    Running standard deviation",
  "              Normalize by n. ",
  "              ",
  "              o(t+(nts-1)/2,x) = std{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runstd1   Running standard deviation (n-1)",
  "              Normalize by (n-1). ",
  "              ",
  "              o(t+(nts-1)/2,x) = std1{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runvar    Running variance",
  "              Normalize by n. ",
  "              ",
  "              o(t+(nts-1)/2,x) = var{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "    runvar1   Running variance (n-1)",
  "              Normalize by (n-1). ",
  "              ",
  "              o(t+(nts-1)/2,x) = var1{i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "",
  "PARAMETER",
  "    nts  INTEGER  Number of timesteps",
  "",
  "ENVIRONMENT",
  "    CDO_TIMESTAT_DATE",
  "        Sets the time stamp in outfile to the \"first\", \"middle\" or "
  "\"last\" contributing timestep of infile.",
};

std::vector<std::string> RunpctlHelp = {
  "NAME",
  "    runpctl - Running percentile values",
  "",
  "SYNOPSIS",
  "    runpctl,p,nts  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes running percentiles over a selected number of "
  "timesteps in infile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "    ",
  "    o(t+(nts-1)/2,x) = pth percentile {i(t,x), i(t+1,x), ..., i(t+nts-1,x)}",
  "",
  "PARAMETER",
  "    p    FLOAT    Percentile number in {0, ..., 100}",
  "    nts  INTEGER  Number of timesteps",
};

std::vector<std::string> TimstatHelp = {
  "NAME",
  "    timmin, timmax, timrange, timsum, timmean, timavg, timstd, timstd1, "
  "timvar, ",
  "    timvar1 - Statistical values over all timesteps",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over all timesteps in infile. "
  "Depending on ",
  "    the chosen operator the minimum, maximum, range, sum, average, variance "
  "or standard deviation of ",
  "    all timesteps read from infile is written to outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    timmin    Time minimum",
  "              o(1,x) = min{i(t',x), t_1<t'<=t_n}",
  "    timmax    Time maximum",
  "              o(1,x) = max{i(t',x), t_1<t'<=t_n}",
  "    timrange  Time range",
  "              o(1,x) = range{i(t',x), t_1<t'<=t_n}",
  "    timsum    Time sum",
  "              o(1,x) = sum{i(t',x), t_1<t'<=t_n}",
  "    timmean   Time mean",
  "              o(1,x) = mean{i(t',x), t_1<t'<=t_n}",
  "    timavg    Time average",
  "              o(1,x) = avg{i(t',x), t_1<t'<=t_n}",
  "    timstd    Time standard deviation",
  "              Normalize by n. ",
  "              ",
  "              o(1,x) = std{i(t',x), t_1<t'<=t_n}",
  "    timstd1   Time standard deviation (n-1)",
  "              Normalize by (n-1). ",
  "              ",
  "              o(1,x) = std1{i(t',x), t_1<t'<=t_n}",
  "    timvar    Time variance",
  "              Normalize by n. ",
  "              ",
  "              o(1,x) = var{i(t',x), t_1<t'<=t_n}",
  "    timvar1   Time variance (n-1)",
  "              Normalize by (n-1). ",
  "              ",
  "              o(1,x) = var1{i(t',x), t_1<t'<=t_n}",
};

std::vector<std::string> TimpctlHelp = {
  "NAME",
  "    timpctl - Percentile values over all timesteps",
  "",
  "SYNOPSIS",
  "    timpctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentiles over all timesteps in infile1. The "
  "algorithm uses ",
  "    histograms with minimum and maximum bounds given in infile2 and "
  "infile3, respectively. ",
  "    The default number of histogram bins is 101. The default can be "
  "overridden by defining the",
  "    environment variable CDO_PCTL_NBINS. The files infile2 and infile3 "
  "should be",
  "    the result of corresponding timmin and timmax operations, respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    ",
  "    o(1,x) = pth percentile {i(t',x), t_1<t'<=t_n}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> HourstatHelp = {
  "NAME",
  "    hourmin, hourmax, hourrange, hoursum, hourmean, houravg, hourstd, "
  "hourstd1, ",
  "    hourvar, hourvar1 - Hourly statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over timesteps of the same "
  "hour.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of timesteps of the same hour is written to "
  "outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    hourmin    Hourly minimum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same hour it is:",
  "               ",
  "               o(t,x) = min{i(t',x), t_1<t'<=t_n}",
  "    hourmax    Hourly maximum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same hour it is:",
  "               ",
  "               o(t,x) = max{i(t',x), t_1<t'<=t_n}",
  "    hourrange  Hourly range",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same hour it is:",
  "               ",
  "               o(t,x) = range{i(t',x), t_1<t'<=t_n}",
  "    hoursum    Hourly sum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same hour it is:",
  "               ",
  "               o(t,x) = sum{i(t',x), t_1<t'<=t_n}",
  "    hourmean   Hourly mean",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same hour it is:",
  "               ",
  "               o(t,x) = mean{i(t',x), t_1<t'<=t_n}",
  "    houravg    Hourly average",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same hour it is:",
  "               ",
  "               o(t,x) = avg{i(t',x), t_1<t'<=t_n}",
  "    hourstd    Hourly standard deviation",
  "               Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same hour it is:",
  "               ",
  "               o(t,x) = std{i(t',x), t_1<t'<=t_n}",
  "    hourstd1   Hourly standard deviation (n-1)",
  "               Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same hour it is:",
  "               ",
  "               o(t,x) = std1{i(t',x), t_1<t'<=t_n}",
  "    hourvar    Hourly variance",
  "               Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same hour it is:",
  "               ",
  "               o(t,x) = var{i(t',x), t_1<t'<=t_n}",
  "    hourvar1   Hourly variance (n-1)",
  "               Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same hour it is:",
  "               ",
  "               o(t,x) = var1{i(t',x), t_1<t'<=t_n}",
};

std::vector<std::string> HourpctlHelp = {
  "NAME",
  "    hourpctl - Hourly percentile values",
  "",
  "SYNOPSIS",
  "    hourpctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentiles over all timesteps of the same hour "
  "in infile1.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and",
  "    infile3, respectively. The default number of histogram bins is 101.",
  "    The default can be overridden by defining the environment variable "
  "CDO_PCTL_NBINS.",
  "    The files infile2 and infile3 should be the result of corresponding "
  "hourmin",
  "    and hourmax operations, respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    For every adjacent sequence t_1, ...,t_n of timesteps of the same hour "
  "it is:",
  "    ",
  "    o(t,x) = pth percentile {i(t',x), t_1<t'<=t_n}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> DaystatHelp = {
  "NAME",
  "    daymin, daymax, dayrange, daysum, daymean, dayavg, daystd, daystd1, "
  "dayvar, ",
  "    dayvar1 - Daily statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over timesteps of the same day.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of timesteps of the same day is written to "
  "outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    daymin    Daily minimum",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same day it is:",
  "              ",
  "              o(t,x) = min{i(t',x), t_1<t'<=t_n}",
  "    daymax    Daily maximum",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same day it is:",
  "              ",
  "              o(t,x) = max{i(t',x), t_1<t'<=t_n}",
  "    dayrange  Daily range",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same day it is:",
  "              ",
  "              o(t,x) = range{i(t',x), t_1<t'<=t_n}",
  "    daysum    Daily sum",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same day it is:",
  "              ",
  "              o(t,x) = sum{i(t',x), t_1<t'<=t_n}",
  "    daymean   Daily mean",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same day it is:",
  "              ",
  "              o(t,x) = mean{i(t',x), t_1<t'<=t_n}",
  "    dayavg    Daily average",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same day it is:",
  "              ",
  "              o(t,x) = avg{i(t',x), t_1<t'<=t_n}",
  "    daystd    Daily standard deviation",
  "              Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same day it is:",
  "              ",
  "              o(t,x) = std{i(t',x), t_1<t'<=t_n}",
  "    daystd1   Daily standard deviation (n-1)",
  "              Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same day it is:",
  "              ",
  "              o(t,x) = std1{i(t',x), t_1<t'<=t_n}",
  "    dayvar    Daily variance",
  "              Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same day it is:",
  "              ",
  "              o(t,x) = var{i(t',x), t_1<t'<=t_n}",
  "    dayvar1   Daily variance (n-1)",
  "              Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same day it is:",
  "              ",
  "              o(t,x) = var1{i(t',x), t_1<t'<=t_n}",
};

std::vector<std::string> DaypctlHelp = {
  "NAME",
  "    daypctl - Daily percentile values",
  "",
  "SYNOPSIS",
  "    daypctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentiles over all timesteps of the same day "
  "in infile1.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and",
  "    infile3, respectively. The default number of histogram bins is 101.",
  "    The default can be overridden by defining the environment variable "
  "CDO_PCTL_NBINS.",
  "    The files infile2 and infile3 should be the result of corresponding "
  "daymin",
  "    and daymax operations, respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    For every adjacent sequence t_1, ...,t_n of timesteps of the same day "
  "it is:",
  "    ",
  "    o(t,x) = pth percentile {i(t',x), t_1<t'<=t_n}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> MonstatHelp = {
  "NAME",
  "    monmin, monmax, monrange, monsum, monmean, monavg, monstd, monstd1, "
  "monvar, ",
  "    monvar1 - Monthly statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over timesteps of the same "
  "month.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of timesteps of the same month is written to "
  "outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    monmin    Monthly minimum",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same month it is:",
  "              ",
  "              o(t,x) = min{i(t',x), t_1<t'<=t_n}",
  "    monmax    Monthly maximum",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same month it is:",
  "              ",
  "              o(t,x) = max{i(t',x), t_1<t'<=t_n}",
  "    monrange  Monthly range",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same month it is:",
  "              ",
  "              o(t,x) = range{i(t',x), t_1<t'<=t_n}",
  "    monsum    Monthly sum",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same month it is:",
  "              ",
  "              o(t,x) = sum{i(t',x), t_1<t'<=t_n}",
  "    monmean   Monthly mean",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same month it is:",
  "              ",
  "              o(t,x) = mean{i(t',x), t_1<t'<=t_n}",
  "    monavg    Monthly average",
  "              For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same month it is:",
  "              ",
  "              o(t,x) = avg{i(t',x), t_1<t'<=t_n}",
  "    monstd    Monthly standard deviation",
  "              Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same month it is:",
  "              ",
  "              o(t,x) = std{i(t',x), t_1 < t' <= t_n}",
  "    monstd1   Monthly standard deviation (n-1)",
  "              Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same month it is:",
  "              ",
  "              o(t,x) = std1{i(t',x), t_1 < t' <= t_n}",
  "    monvar    Monthly variance",
  "              Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same month it is:",
  "              ",
  "              o(t,x) = var{i(t',x), t_1 < t' <= t_n}",
  "    monvar1   Monthly variance (n-1)",
  "              Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same month it is:",
  "              ",
  "              o(t,x) = var1{i(t',x), t_1 < t' <= t_n}",
};

std::vector<std::string> MonpctlHelp = {
  "NAME",
  "    monpctl - Monthly percentile values",
  "",
  "SYNOPSIS",
  "    monpctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentiles over all timesteps of the same month "
  "in infile1.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and",
  "    infile3, respectively. The default number of histogram bins is 101.",
  "    The default can be overridden by defining the environment variable "
  "CDO_PCTL_NBINS.",
  "    The files infile2 and infile3 should be the result of corresponding "
  "monmin",
  "    and monmax operations, respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    For every adjacent sequence t_1, ...,t_n of timesteps of the same month "
  "it is:",
  "    ",
  "    o(t,x) = pth percentile {i(t',x), t_1<t'<=t_n}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> YearmonstatHelp = {
  "NAME",
  "    yearmonmean - Yearly mean from monthly data",
  "",
  "SYNOPSIS",
  "    yearmonmean  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator computes the yearly mean of a monthly time series.",
  "    Each month is weighted with the number of days per month. ",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "    ",
  "    For every adjacent sequence t_1, ...,t_n of timesteps of the same year "
  "it is:",
  "    ",
  "    o(t,x) = mean{i(t',x), t_1<t'<=t_n}",
  "",
  "ENVIRONMENT",
  "    CDO_TIMESTAT_DATE",
  "        Sets the date information in outfile to the \"first\", \"middle\" "
  "or \"last\" contributing timestep of infile.",
};

std::vector<std::string> YearstatHelp = {
  "NAME",
  "    yearmin, yearmax, yearrange, yearsum, yearmean, yearavg, yearstd, "
  "yearstd1, ",
  "    yearvar, yearvar1 - Yearly statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over timesteps of the same "
  "year.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of timesteps of the same year is written to "
  "outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "",
  "OPERATORS",
  "    yearmin    Yearly minimum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same year it is:",
  "               ",
  "               o(t,x) = min{i(t',x), t_1<t'<=t_n}",
  "    yearmax    Yearly maximum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same year it is:",
  "               ",
  "               o(t,x) = max{i(t',x), t_1<t'<=t_n}",
  "    yearrange  Yearly range",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same year it is:",
  "               ",
  "               o(t,x) = range{i(t',x), t_1<t'<=t_n}",
  "    yearsum    Yearly sum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same year it is:",
  "               ",
  "               o(t,x) = sum{i(t',x), t_1<t'<=t_n}",
  "    yearmean   Yearly mean",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same year it is:",
  "               ",
  "               o(t,x) = mean{i(t',x), t_1<t'<=t_n}",
  "    yearavg    Yearly average",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same year it is:",
  "               ",
  "               o(t,x) = avg{i(t',x), t_1<t'<=t_n}",
  "    yearstd    Yearly standard deviation",
  "               Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same year it is:",
  "               ",
  "               o(t,x) = std{i(t',x), t_1 < t' <= t_n}",
  "    yearstd1   Yearly standard deviation (n-1)",
  "               Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same year it is:",
  "               ",
  "               o(t,x) = std1{i(t',x), t_1 < t' <= t_n}",
  "    yearvar    Yearly variance",
  "               Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same year it is:",
  "               ",
  "               o(t,x) = var{i(t',x), t_1 < t' <= t_n}",
  "    yearvar1   Yearly variance (n-1)",
  "               Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same year it is:",
  "               ",
  "               o(t,x) = var1{i(t',x), t_1 < t' <= t_n}",
  "",
  "NOTE",
  "    The operators yearmean and yearavg compute only arithmetical means!",
};

std::vector<std::string> YearpctlHelp = {
  "NAME",
  "    yearpctl - Yearly percentile values",
  "",
  "SYNOPSIS",
  "    yearpctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentiles over all timesteps of the same year "
  "in infile1.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and",
  "    infile3, respectively. The default number of histogram bins is 101. The "
  "default can be",
  "    overridden by defining the environment variable CDO_PCTL_NBINS. The "
  "files infile2 and",
  "    infile3 should be the result of corresponding yearmin and yearmax "
  "operations, respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    For every adjacent sequence t_1, ...,t_n of timesteps of the same year "
  "it is:",
  "    ",
  "    o(t,x) = pth percentile {i(t',x), t_1<t'<=t_n}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> SeasstatHelp = {
  "NAME",
  "    seasmin, seasmax, seasrange, seassum, seasmean, seasavg, seasstd, "
  "seasstd1, ",
  "    seasvar, seasvar1 - Seasonal statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values over timesteps of the same "
  "season.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of timesteps of the same season is written to "
  "outfile.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile.",
  "    Be careful about the first and the last output timestep, they may be "
  "incorrect values ",
  "    if the seasons have incomplete timesteps.",
  "",
  "OPERATORS",
  "    seasmin    Seasonal minimum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same season it is:",
  "               ",
  "               o(t,x) = min{i(t',x), t1 < t' <= tn}",
  "    seasmax    Seasonal maximum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same season it is:",
  "               ",
  "               o(t,x) = max{i(t',x), t1 < t' <= tn}",
  "    seasrange  Seasonal range",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same season it is:",
  "               ",
  "               o(t,x) = range{i(t',x), t1 < t' <= tn}",
  "    seassum    Seasonal sum",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same season it is:",
  "               ",
  "               o(t,x) = sum{i(t',x), t1 < t' <= tn}",
  "    seasmean   Seasonal mean",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same season it is:",
  "               ",
  "               o(t,x) = mean{i(t',x), t1 < t' <= tn}",
  "    seasavg    Seasonal average",
  "               For every adjacent sequence t_1, ...,t_n of timesteps of the "
  "same season it is:",
  "               ",
  "               o(t,x) = avg{i(t',x), t1 < t' <= tn}",
  "    seasstd    Seasonal standard deviation",
  "               Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same season it is:",
  "               ",
  "               o(t,x) = std{i(t',x), t1 < t' <= tn}",
  "    seasstd1   Seasonal standard deviation (n-1)",
  "               Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same season it is:",
  "               ",
  "               o(t,x) = std1{i(t',x), t1 < t' <= tn}",
  "    seasvar    Seasonal variance",
  "               Normalize by n. For every adjacent sequence t_1, ...,t_n of "
  "timesteps of the same season it is:",
  "               ",
  "               o(t,x) = var{i(t',x), t1 < t' <= tn}",
  "    seasvar1   Seasonal variance (n-1)",
  "               Normalize by (n-1). For every adjacent sequence t_1, ...,t_n "
  "of timesteps of the same season it is:",
  "               ",
  "               o(t,x) = var1{i(t',x), t1 < t' <= tn}",
};

std::vector<std::string> SeaspctlHelp = {
  "NAME",
  "    seaspctl - Seasonal percentile values",
  "",
  "SYNOPSIS",
  "    seaspctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator computes percentiles over all timesteps in infile1 of the "
  "same season.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and infile3,",
  "    respectively. The default number of histogram bins is 101. The default "
  "can be overridden",
  "    by defining the environment variable CDO_PCTL_NBINS. The files infile2 "
  "and infile3",
  "    should be the result of corresponding seasmin and seasmax operations, "
  "respectively.",
  "    The time of outfile is determined by the time in the middle of all "
  "contributing timesteps of infile1.",
  "    Be careful about the first and the last output timestep, they may be "
  "incorrect values ",
  "    if the seasons have incomplete timesteps.",
  "    For every adjacent sequence t_1, ...,t_n of timesteps of the same "
  "season it is:",
  "    ",
  "    o(t,x) = pth percentile {i(t',x), t1 < t' <= tn}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> YhourstatHelp = {
  "NAME",
  "    yhourmin, yhourmax, yhourrange, yhoursum, yhourmean, yhouravg, "
  "yhourstd, ",
  "    yhourstd1, yhourvar, yhourvar1 - Multi-year hourly statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values of each hour and day of year.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of each hour and day of year in infile is written "
  "to outfile.",
  "    The date information in an output field is the date of the last "
  "contributing input field.",
  "",
  "OPERATORS",
  "    yhourmin    Multi-year hourly minimum",
  "                o(0001,x) = min{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = min{i(t,x), day(i(t)) = 8784}",
  "    yhourmax    Multi-year hourly maximum",
  "                o(0001,x) = max{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = max{i(t,x), day(i(t)) = 8784}",
  "    yhourrange  Multi-year hourly range",
  "                o(0001,x) = range{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = range{i(t,x), day(i(t)) = 8784}",
  "    yhoursum    Multi-year hourly sum",
  "                o(0001,x) = sum{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = sum{i(t,x), day(i(t)) = 8784}",
  "    yhourmean   Multi-year hourly mean",
  "                o(0001,x) = mean{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = mean{i(t,x), day(i(t)) = 8784}",
  "    yhouravg    Multi-year hourly average",
  "                o(0001,x) = avg{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = avg{i(t,x), day(i(t)) = 8784}",
  "    yhourstd    Multi-year hourly standard deviation",
  "                Normalize by n. ",
  "                ",
  "                o(0001,x) = std{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = std{i(t,x), day(i(t)) = 8784}",
  "    yhourstd1   Multi-year hourly standard deviation (n-1)",
  "                Normalize by (n-1). ",
  "                ",
  "                o(0001,x) = std1{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = std1{i(t,x), day(i(t)) = 8784}",
  "    yhourvar    Multi-year hourly variance",
  "                Normalize by n. ",
  "                ",
  "                o(0001,x) = var{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = var{i(t,x), day(i(t)) = 8784}",
  "    yhourvar1   Multi-year hourly variance (n-1)",
  "                Normalize by (n-1). ",
  "                ",
  "                o(0001,x) = var1{i(t,x), day(i(t)) = 0001}",
  "                                 ...",
  "                o(8784,x) = var1{i(t,x), day(i(t)) = 8784}",
};

std::vector<std::string> YdaystatHelp = {
  "NAME",
  "    ydaymin, ydaymax, ydayrange, ydaysum, ydaymean, ydayavg, ydaystd, "
  "ydaystd1, ",
  "    ydayvar, ydayvar1 - Multi-year daily statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values of each day of year.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of each day of year in infile is written to "
  "outfile.",
  "    The date information in an output field is the date of the last "
  "contributing input field.",
  "",
  "OPERATORS",
  "    ydaymin    Multi-year daily minimum",
  "               o(001,x) = min{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = min{i(t,x), day(i(t)) = 366}",
  "    ydaymax    Multi-year daily maximum",
  "               o(001,x) = max{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = max{i(t,x), day(i(t)) = 366}",
  "    ydayrange  Multi-year daily range",
  "               o(001,x) = range{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = range{i(t,x), day(i(t)) = 366}",
  "    ydaysum    Multi-year daily sum",
  "               o(001,x) = sum{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = sum{i(t,x), day(i(t)) = 366}",
  "    ydaymean   Multi-year daily mean",
  "               o(001,x) = mean{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = mean{i(t,x), day(i(t)) = 366}",
  "    ydayavg    Multi-year daily average",
  "               o(001,x) = avg{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = avg{i(t,x), day(i(t)) = 366}",
  "    ydaystd    Multi-year daily standard deviation",
  "               Normalize by n. ",
  "               ",
  "               o(001,x) = std{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = std{i(t,x), day(i(t)) = 366}",
  "    ydaystd1   Multi-year daily standard deviation (n-1)",
  "               Normalize by (n-1). ",
  "               ",
  "               o(001,x) = std1{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = std1{i(t,x), day(i(t)) = 366}",
  "    ydayvar    Multi-year daily variance",
  "               Normalize by n. ",
  "               ",
  "               o(001,x) = var{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = var{i(t,x), day(i(t)) = 366}",
  "    ydayvar1   Multi-year daily variance (n-1)",
  "               Normalize by (n-1). ",
  "               ",
  "               o(001,x) = var1{i(t,x), day(i(t)) = 001}",
  "                                ...",
  "               o(366,x) = var1{i(t,x), day(i(t)) = 366}",
};

std::vector<std::string> YdaypctlHelp = {
  "NAME",
  "    ydaypctl - Multi-year daily percentile values",
  "",
  "SYNOPSIS",
  "    ydaypctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator writes a certain percentile of each day of year in "
  "infile1 to outfile.",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and",
  "    infile3, respectively. The default number of histogram bins is 101. The "
  "default can be",
  "    overridden by setting the environment variable CDO_PCTL_NBINS to a "
  "different value.",
  "    The files infile2 and infile3 should be the result of corresponding "
  "ydaymin and",
  "    ydaymax operations, respectively.",
  "    The date information in an output field is the date of the last "
  "contributing input field.",
  "    ",
  "    o(001,x) = pth percentile {i(t,x), day(i(t)) = 001}",
  "                     ...",
  "    o(366,x) = pth percentile {i(t,x), day(i(t)) = 366}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> YmonstatHelp = {
  "NAME",
  "    ymonmin, ymonmax, ymonrange, ymonsum, ymonmean, ymonavg, ymonstd, "
  "ymonstd1, ",
  "    ymonvar, ymonvar1 - Multi-year monthly statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values of each month of year.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of each month of year in infile is written to "
  "outfile.",
  "    The date information in an output field is the date of the last "
  "contributing input field.",
  "",
  "OPERATORS",
  "    ymonmin    Multi-year monthly minimum",
  "               o(01,x) = min{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = min{i(t,x), month(i(t)) = 12}",
  "    ymonmax    Multi-year monthly maximum",
  "               o(01,x) = max{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = max{i(t,x), month(i(t)) = 12}",
  "    ymonrange  Multi-year monthly range",
  "               o(01,x) = range{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = range{i(t,x), month(i(t)) = 12}",
  "    ymonsum    Multi-year monthly sum",
  "               o(01,x) = sum{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = sum{i(t,x), month(i(t)) = 12}",
  "    ymonmean   Multi-year monthly mean",
  "               o(01,x) = mean{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = mean{i(t,x), month(i(t)) = 12}",
  "    ymonavg    Multi-year monthly average",
  "               o(01,x) = avg{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = avg{i(t,x), month(i(t)) = 12}",
  "    ymonstd    Multi-year monthly standard deviation",
  "               Normalize by n. ",
  "               ",
  "               o(01,x) = std{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = std{i(t,x), month(i(t)) = 12}",
  "    ymonstd1   Multi-year monthly standard deviation (n-1)",
  "               Normalize by (n-1). ",
  "               ",
  "               o(01,x) = std1{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = std1{i(t,x), month(i(t)) = 12}",
  "    ymonvar    Multi-year monthly variance",
  "               Normalize by n. ",
  "               ",
  "               o(01,x) = var{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = var{i(t,x), month(i(t)) = 12}",
  "    ymonvar1   Multi-year monthly variance (n-1)",
  "               Normalize by (n-1). ",
  "               ",
  "               o(01,x) = var1{i(t,x), month(i(t)) = 01}",
  "                                ...",
  "               o(12,x) = var1{i(t,x), month(i(t)) = 12}",
};

std::vector<std::string> YmonpctlHelp = {
  "NAME",
  "    ymonpctl - Multi-year monthly percentile values",
  "",
  "SYNOPSIS",
  "    ymonpctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator writes a certain percentile of each month of year in "
  "infile1 to outfile.",
  "    The algorithm uses histograms with minimum and maximum bounds given in",
  "    infile2 and infile3, respectively. The default number of",
  "    histogram bins is 101. The default can be overridden by setting the",
  "    environment variable CDO_PCTL_NBINS to a different value. The files",
  "    infile2 and infile3 should be the result of corresponding",
  "    ymonmin and ymonmax operations, respectively.",
  "    The date information in an output field is the date of the last",
  "    contributing input field.",
  "    ",
  "    o(01,x) = pth percentile {i(t,x), month(i(t)) = 01}",
  "                     ...",
  "    o(12,x) = pth percentile {i(t,x), month(i(t)) = 12}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> YseasstatHelp = {
  "NAME",
  "    yseasmin, yseasmax, yseasrange, yseassum, yseasmean, yseasavg, "
  "yseasstd, ",
  "    yseasstd1, yseasvar, yseasvar1 - Multi-year seasonal statistical values",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module computes statistical values of each season.",
  "    Depending on the chosen operator the minimum, maximum, range, sum, "
  "average, variance",
  "    or standard deviation of each season in infile is written to outfile.",
  "    The date information in an output field is the date of the last "
  "contributing input field.",
  "",
  "OPERATORS",
  "    yseasmin    Multi-year seasonal minimum",
  "                o(1,x) = min{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = min{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = min{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = min{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasmax    Multi-year seasonal maximum",
  "                o(1,x) = max{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = max{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = max{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = max{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasrange  Multi-year seasonal range",
  "                o(1,x) = range{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = range{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = range{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = range{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseassum    Multi-year seasonal sum",
  "                o(1,x) = sum{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = sum{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = sum{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = sum{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasmean   Multi-year seasonal mean",
  "                o(1,x) = mean{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = mean{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = mean{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = mean{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasavg    Multi-year seasonal average",
  "                o(1,x) = avg{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = avg{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = avg{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = avg{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasstd    Multi-year seasonal standard deviation",
  "                o(1,x) = std{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = std{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = std{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = std{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasstd1   Multi-year seasonal standard deviation (n-1)",
  "                o(1,x) = std1{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = std1{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = std1{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = std1{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasvar    Multi-year seasonal variance",
  "                o(1,x) = var{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = var{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = var{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = var{i(t,x), month(i(t)) = 09, 10, 11}",
  "    yseasvar1   Multi-year seasonal variance (n-1)",
  "                o(1,x) = var1{i(t,x), month(i(t)) = 12, 01, 02}",
  "                o(2,x) = var1{i(t,x), month(i(t)) = 03, 04, 05}",
  "                o(3,x) = var1{i(t,x), month(i(t)) = 06, 07, 08}",
  "                o(4,x) = var1{i(t,x), month(i(t)) = 09, 10, 11}",
};

std::vector<std::string> YseaspctlHelp = {
  "NAME",
  "    yseaspctl - Multi-year seasonal percentile values",
  "",
  "SYNOPSIS",
  "    yseaspctl,p  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator writes a certain percentile of each season in infile1 to "
  "outfile.",
  "    The algorithm uses histograms with minimum and maximum bounds given in",
  "    infile2 and infile3, respectively. The default number of",
  "    histogram bins is 101. The default can be overridden by setting the",
  "    environment variable CDO_PCTL_NBINS to a different value. The files",
  "    infile2 and infile3 should be the result of corresponding",
  "    yseasmin and yseasmax operations, respectively.",
  "    The date information in an output field is the date of the last",
  "    contributing input field.",
  "    ",
  "    o(1,x) = pth percentile {i(t,x), month(i(t)) = 12, 01, 02}",
  "    o(2,x) = pth percentile {i(t,x), month(i(t)) = 03, 04, 05}",
  "    o(3,x) = pth percentile {i(t,x), month(i(t)) = 06, 07, 08}",
  "    o(4,x) = pth percentile {i(t,x), month(i(t)) = 09, 10, 11}",
  "",
  "PARAMETER",
  "    p  FLOAT  Percentile number in {0, ..., 100}",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> YdrunstatHelp = {
  "NAME",
  "    ydrunmin, ydrunmax, ydrunsum, ydrunmean, ydrunavg, ydrunstd, "
  "ydrunstd1, ",
  "    ydrunvar, ydrunvar1 - Multi-year daily running statistical values",
  "",
  "SYNOPSIS",
  "    <operator>,nts  infile outfile",
  "",
  "DESCRIPTION",
  "    This module writes running statistical values for each day of year in "
  "infile to outfile.",
  "    Depending on the chosen operator, the minimum, maximum, sum, average, "
  "variance or standard deviation ",
  "    of all timesteps in running windows of which the medium timestep "
  "corresponds to a certain day of",
  "    year is computed. The date information in an output field is the date "
  "of the timestep in the middle ",
  "    of the last contributing running window.",
  "    Note that the operator have to be applied to a continuous time series "
  "of daily measurements in order ",
  "    to yield physically meaningful results. Also note that the output time "
  "series begins (nts-1)/2 timesteps",
  "    after the first timestep of the input time series and ends (nts-1)/2 "
  "timesteps before the last one.",
  "    For input data which are complete but not continuous, such as time "
  "series of daily measurements for ",
  "    the same month or season within different years, the operator yields "
  "physically meaningful results ",
  "    only if the input time series does include the (nts-1)/2 days before "
  "and after each period of interest.",
  "",
  "OPERATORS",
  "    ydrunmin   Multi-year daily running minimum",
  "               o(001,x) = min{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = min{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "    ydrunmax   Multi-year daily running maximum",
  "               o(001,x) = max{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = max{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "    ydrunsum   Multi-year daily running sum",
  "               o(001,x) = sum{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = sum{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "    ydrunmean  Multi-year daily running mean",
  "               o(001,x) = mean{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = mean{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "    ydrunavg   Multi-year daily running average",
  "               o(001,x) = avg{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = avg{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "    ydrunstd   Multi-year daily running standard deviation",
  "               Normalize by n. ",
  "               ",
  "               o(001,x) = std{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = std{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[i(t+(nts-1)/2)] = 366}",
  "    ydrunstd1  Multi-year daily running standard deviation (n-1)",
  "               Normalize by (n-1). ",
  "               ",
  "               o(001,x) = std1{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = std1{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[i(t+(nts-1)/2)] = 366}",
  "    ydrunvar   Multi-year daily running variance",
  "               Normalize by n. ",
  "               ",
  "               o(001,x) = var{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = var{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "    ydrunvar1  Multi-year daily running variance (n-1)",
  "               Normalize by (n-1). ",
  "               ",
  "               o(001,x) = var1{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                                ...",
  "               o(366,x) = var1{i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "",
  "PARAMETER",
  "    nts  INTEGER  Number of timesteps",
};

std::vector<std::string> YdrunpctlHelp = {
  "NAME",
  "    ydrunpctl - Multi-year daily running percentile values",
  "",
  "SYNOPSIS",
  "    ydrunpctl,p,nts  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator writes running percentile values for each day of year in "
  "infile1 to outfile. ",
  "    A certain percentile is computed for all timesteps in running windows "
  "of which the medium ",
  "    timestep corresponds to a certain day of year. ",
  "    The algorithm uses histograms with minimum and maximum bounds given in "
  "infile2 and infile3,",
  "    respectively. The default number of histogram bins is 101. The default "
  "can be overridden",
  "    by setting the environment variable CDO_PCTL_NBINS to a different "
  "value. The files infile2 ",
  "    and infile3 should be the result of corresponding ydrunmin and ydrunmax "
  "operations, respectively.",
  "    The date information in an output field is the date of the timestep in "
  "the middle of the last ",
  "    contributing running window.",
  "    Note that the operator have to be applied to a continuous time series "
  "of daily measurements ",
  "    in order to yield physically meaningful results. Also note that the "
  "output time series begins",
  "    (nts-1)/2 timesteps after the first timestep of the input time series "
  "and ends (nts-1)/2 ",
  "    timesteps before the last.",
  "    For input data which are complete but not continuous, such as time "
  "series of daily measurements ",
  "    for the same month or season within different years, the operator only "
  "yields physically meaningful ",
  "    results if the input time series does include the (nts-1)/2 days before "
  "and after each period ",
  "    of interest.",
  "    ",
  "    o(001,x) = pth percentile {i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 001}",
  "                     ...",
  "    o(366,x) = pth percentile {i(t,x), i(t+1,x), ..., i(t+nts-1,x); "
  "day[(i(t+(nts-1)/2)] = 366}",
  "",
  "PARAMETER",
  "    p    FLOAT    Percentile number in {0, ..., 100}",
  "    nts  INTEGER  Number of timesteps",
  "",
  "ENVIRONMENT",
  "    CDO_PCTL_NBINS",
  "        Sets the number of histogram bins. The default number is 101.",
};

std::vector<std::string> FldcorHelp = {
  "NAME",
  "    fldcor - Correlation in grid space",
  "",
  "SYNOPSIS",
  "    fldcor  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    The correlation coefficient is a quantity that gives the quality of a "
  "least ",
  "    squares fitting to the original data. This operator correlates all "
  "gridpoints",
  "    of two fields for each timestep. With",
  "    ",
  "    S(t) = {x, i_1(t,x) != missval and i_2(t,x) != missval}",
  "    ",
  "    it is",
  "    ",
  "    o(t,1) = Cor{(i_1(t,x), i_2(t,x)), x_1 < x <= x_n}",
  "    ",
  "    where w(x) are the area weights obtained by the input streams.",
  "    For every timestep t only those field elements x belong to the sample,",
  "    which have i_1(t,x) != missval and i_2(t,x) != missval.",
};

std::vector<std::string> TimcorHelp = {
  "NAME",
  "    timcor - Correlation over time",
  "",
  "SYNOPSIS",
  "    timcor  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    The correlation coefficient is a quantity that gives the quality of a "
  "least ",
  "    squares fitting to the original data. This operator correlates each "
  "gridpoint",
  "    of two fields over all timesteps. With",
  "    ",
  "    S(x) = {t, i_1(t,x) != missval and i_2(t,x) != missval}",
  "    ",
  "    it is",
  "    ",
  "    o(1,x) = Cor{(i_1(t,x), i_2(t,x)), t_1 < t <= t_n}",
  "    ",
  "    For every gridpoint x only those timesteps t belong to the sample,",
  "    which have i_1(t,x) != missval and i_2(t,x) != missval.",
};

std::vector<std::string> FldcovarHelp = {
  "NAME",
  "    fldcovar - Covariance in grid space",
  "",
  "SYNOPSIS",
  "    fldcovar  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This operator calculates the covariance of two fields over all "
  "gridpoints",
  "    for each timestep. With",
  "    ",
  "    S(t) = {x, i_1(t,x) != missval and i_2(t,x) != missval}",
  "    ",
  "    it is",
  "    ",
  "    o(t,1) = Covar{(i_1(t,x), i_2(t,x)), x_1 < x <= x_n}",
  "    ",
  "    where w(x) are the area weights obtained by the input streams.",
  "    For every timestep t only those field elements x belong to the sample,",
  "    which have i_1(t,x) != missval and i_2(t,x) != missval.",
};

std::vector<std::string> TimcovarHelp = {
  "NAME",
  "    timcovar - Covariance over time",
  "",
  "SYNOPSIS",
  "    timcovar  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This operator calculates the covariance of two fields at each gridpoint",
  "    over all timesteps. With",
  "    ",
  "    S(x) = {t, i_1(t,x) != missval and i_2(t,x) != missval}",
  "    ",
  "    it is",
  "    ",
  "    o(1,x) = Covar{(i_1(t,x), i_2(t,x)), t_1 < t <= t_n}",
  "    ",
  "    For every gridpoint x only those timesteps t belong to the sample,",
  "    which have i_1(t,x) != missval and i_2(t,x) != missval.",
};

std::vector<std::string> RegresHelp = {
  "NAME",
  "    regres - Regression",
  "",
  "SYNOPSIS",
  "    regres  infile outfile",
  "",
  "DESCRIPTION",
  "    The values of the input file infile are assumed to be distributed as",
  "    N(a+b*t,S^2) with unknown a, b and S^2. This operator estimates the",
  "    parameter b. For every field element x only those timesteps ",
  "    t belong to the sample S(x), which have i(t,x) NE miss.",
};

std::vector<std::string> DetrendHelp = {
  "NAME",
  "    detrend - Detrend time series",
  "",
  "SYNOPSIS",
  "    detrend  infile outfile",
  "",
  "DESCRIPTION",
  "    Every time series in infile is linearly detrended. For every field "
  "element x ",
  "    only those timesteps t belong to the sample S(x), which have i(t,x) NE "
  "miss.",
  "",
  "NOTE",
  "    This operator has to keep the fields of all timesteps concurrently in "
  "the memory.",
  "    If not enough memory is available use the operators trend and subtrend.",
};

std::vector<std::string> TrendHelp = {
  "NAME",
  "    trend - Trend of time series",
  "",
  "SYNOPSIS",
  "    trend  infile outfile1 outfile2",
  "",
  "DESCRIPTION",
  "    The values of the input file infile are assumed to be distributed as",
  "    N(a+b*t,S^2) with unknown a, b and S^2. This operator estimates the",
  "    parameter a and b. For every field element x only those timesteps ",
  "    t belong to the sample S(x), which have i(t,x) NE miss.",
  "    Thus the estimation for a is stored in outfile1 and that for b is "
  "stored ",
  "    in outfile2. To subtract the trend from the data see operator subtrend.",
};

std::vector<std::string> SubtrendHelp = {
  "NAME",
  "    subtrend - Subtract a trend",
  "",
  "SYNOPSIS",
  "    subtrend  infile1 infile2 infile3 outfile",
  "",
  "DESCRIPTION",
  "    This operator is for subtracting a trend computed by the operator "
  "trend.",
  "    It is",
  "    ",
  "    o(t,x) = i_1(t,x) - (i_2(1,x) + i_3(1,x)*t)",
  "    where t is the timesteps.",
};

std::vector<std::string> EOFsHelp = {
  "NAME",
  "    eof, eoftime, eofspatial, eof3d - Empirical Orthogonal Functions",
  "",
  "SYNOPSIS",
  "    <operator>,neof  infile outfile1 outfile2",
  "",
  "DESCRIPTION",
  "    This module calculates empirical orthogonal functions of the data in "
  "infile ",
  "    as the eigen values of the scatter matrix (covariance matrix) S of the "
  "data",
  "    sample z(t). A more detailed description can be found above.",
  "    ",
  "    Please note, that the input data are assumed to be anomalies.",
  "    ",
  "    If operator eof is chosen, the EOFs are computed in either time or "
  "spatial",
  "    space, whichever is the fastest. If the user already knows, which "
  "computation",
  "    is faster, the module can be forced to perform a computation in time- "
  "or gridspace",
  "    by using the operators eoftime or eofspatial, respectively. This can "
  "enhance ",
  "    performance, especially for very long time series, where the number of "
  "timesteps",
  "    is larger than the number of grid-points. Data in infile are assumed to "
  "be anomalies.",
  "    If they are not, the behavior of this module is not well defined. ",
  "    After execution outfile1 will contain all eigen-values and outfile2 the",
  "    eigenvectors e_j. All EOFs and eigen-values are computed. However, only "
  "the first ",
  "    neof EOFs are written to outfile2. Nonetheless, outfile1 contains all "
  "eigen-values. ",
  "    ",
  "    Missing values are not fully supported. Support is only checked for "
  "non-changing",
  "    masks of missing values in time. Although there still will be results, "
  "they are",
  "    not trustworthy, and a warning will occur. In the latter case we "
  "suggest to ",
  "    replace missing values by 0 in infile. ",
  "",
  "OPERATORS",
  "    eof         Calculate EOFs in spatial or time space",
  "    eoftime     Calculate EOFs in time space",
  "    eofspatial  Calculate EOFs in spatial space",
  "    eof3d       Calculate 3-Dimensional EOFs in time space",
  "",
  "PARAMETER",
  "    neof  INTEGER  Number of eigen functions",
  "",
  "ENVIRONMENT",
  "    CDO_SVD_MODE   ",
  "        Is used to choose the algorithm for eigenvalue calculation. Options "
  "are 'jacobi' for ",
  "        a one-sided parallel jacobi-algorithm (only executed in parallel if "
  "-P flag is set)",
  "        and  'danielson_lanczos' for a non-parallel d/l algorithm. The "
  "default setting is 'jacobi'.",
  "    CDO_WEIGHT_MODE",
  "        It is used to set the weight mode. The default is 'off'. Set it to "
  "'on' for a weighted version.",
  "    MAX_JACOBI_ITER",
  "        Is the maximum integer number of annihilation sweeps that is "
  "executed if the ",
  "        jacobi-algorithm is used to compute the eigen values. The default "
  "value is 12.",
  "    FNORM_PRECISION",
  "        Is the Frobenius norm of the matrix consisting of an annihilation "
  "pair",
  "        of eigenvectors that is used to determine if the eigenvectors have "
  "reached ",
  "        a sufficient level of convergence. If all annihilation-pairs of "
  "vectors have ",
  "        a norm below this value, the computation is considered to have "
  "converged ",
  "        properly. Otherwise, a warning will occur. The default value 1e-12.",
};

std::vector<std::string> EofcoeffHelp = {
  "NAME",
  "    eofcoeff - Principal coefficients of EOFs",
  "",
  "SYNOPSIS",
  "    eofcoeff  infile1 infile2 obase",
  "",
  "DESCRIPTION",
  "    This module calculates the time series of the principal coefficients "
  "for given EOF",
  "    (empirical orthogonal functions) and data. Time steps in infile1 are "
  "assumed to be the EOFs,",
  "    time steps in infile2 are assumed to be the time series.",
  "    Note, that this operator calculates a non weighted dot product of the "
  "fields in infile1 and infile2.",
  "    For consistency set the environment variable CDO_WEIGHT_MODE=off when "
  "using eof or eof3d.",
  "    ",
  "    There will be a separate file containing a time series of principal "
  "coefficients",
  "    with time information from infile2 for each EOF in infile1. Output "
  "files will be",
  "    numbered as <obase><neof><suffix> where neof+1 is the number of the EOF "
  "(timestep)",
  "    in infile1 and suffix is the filename extension derived from the file "
  "format. ",
  "",
  "ENVIRONMENT",
  "    CDO_FILE_SUFFIX",
  "        Set the default file suffix. This suffix will be added to the "
  "output file ",
  "        names instead of the filename extension derived from the file "
  "format. ",
  "        Set this variable to NULL to disable the adding of a file suffix.",
};

std::vector<std::string> RemapbilHelp = {
  "NAME",
  "    remapbil, genbil - Bilinear interpolation",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a bilinear remapping of fields "
  "between grids in spherical coordinates.",
  "    The interpolation is based on an adapted SCRIP library version. ",
  "    For a detailed description of the interpolation method see SCRIP.",
  "    This interpolation method only works on quadrilateral curvilinear "
  "source grids.",
  "",
  "OPERATORS",
  "    remapbil  Bilinear interpolation",
  "              Performs a bilinear interpolation on all input fields.",
  "    genbil    Generate bilinear interpolation weights",
  "              Generates bilinear interpolation weights for the first input "
  "field and writes the",
  "              result to a file. The format of this file is NetCDF following "
  "the SCRIP convention.",
  "              Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    REMAP_EXTRAPOLATE",
  "        This variable is used to switch the extrapolation feature 'on' or "
  "'off'.",
  "        By default the extrapolation is enabled for circular grids.",
};

std::vector<std::string> RemapbicHelp = {
  "NAME",
  "    remapbic, genbic - Bicubic interpolation",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a bicubic remapping of fields "
  "between grids in spherical coordinates.",
  "    The interpolation is based on an adapted SCRIP library version. ",
  "    For a detailed description of the interpolation method see SCRIP.",
  "    This interpolation method only works on quadrilateral curvilinear "
  "source grids.",
  "",
  "OPERATORS",
  "    remapbic  Bicubic interpolation",
  "              Performs a bicubic interpolation on all input fields.",
  "    genbic    Generate bicubic interpolation weights",
  "              Generates bicubic interpolation weights for the first input "
  "field and writes the",
  "              result to a file. The format of this file is NetCDF following "
  "the SCRIP convention.",
  "              Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    REMAP_EXTRAPOLATE",
  "        This variable is used to switch the extrapolation feature 'on' or "
  "'off'.",
  "        By default the extrapolation is enabled for circular grids.",
};

std::vector<std::string> RemapnnHelp = {
  "NAME",
  "    remapnn, gennn - Nearest neighbor remapping",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a nearest neighbor remapping of "
  "fields between grids",
  "    in spherical coordinates.",
  "",
  "OPERATORS",
  "    remapnn  Nearest neighbor remapping",
  "             Performs a nearest neighbor remapping on all input fields.",
  "    gennn    Generate nearest neighbor remap weights",
  "             Generates nearest neighbor remapping weights for the first "
  "input field and writes the result to a file.",
  "             The format of this file is NetCDF following the SCRIP "
  "convention.",
  "             Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    REMAP_EXTRAPOLATE    ",
  "        This variable is used to switch the extrapolation feature 'on' or "
  "'off'.",
  "        By default the extrapolation is enabled for this remapping method.",
  "    CDO_GRIDSEARCH_RADIUS",
  "        Grid search radius in degree, default 180 degree.",
};

std::vector<std::string> RemapdisHelp = {
  "NAME",
  "    remapdis, gendis - Distance-weighted average remapping",
  "",
  "SYNOPSIS",
  "    remapdis,grid[,neighbors]  infile outfile",
  "    gendis,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a distance-weighted average "
  "remapping of the four",
  "    nearest neighbor values of fields between grids in spherical "
  "coordinates.",
  "    The interpolation is based on an adapted SCRIP library version. ",
  "    For a detailed description of the interpolation method see SCRIP.",
  "",
  "OPERATORS",
  "    remapdis  Distance-weighted average remapping",
  "              Performs a distance-weighted average remapping of the nearest "
  "neighbors value on all input fields.",
  "              The default number of nearest neighbors is 4.",
  "    gendis    Generate distance-weighted average remap weights",
  "              Generates distance-weighted average remapping weights of the "
  "four nearest neighbor",
  "              values for the first input field and writes the result to a "
  "file.",
  "              The format of this file is NetCDF following the SCRIP "
  "convention.",
  "              Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid       STRING   Target grid description file or name",
  "    neighbors  INTEGER  Number of nearest neighbors",
  "",
  "ENVIRONMENT",
  "    REMAP_EXTRAPOLATE    ",
  "        This variable is used to switch the extrapolation feature 'on' or "
  "'off'.",
  "        By default the extrapolation is enabled for this remapping method.",
  "    CDO_GRIDSEARCH_RADIUS",
  "        Grid search radius in degree, default 180 degree.",
};

std::vector<std::string> RemapyconHelp = {
  "NAME",
  "    remapycon, genycon - First order conservative remapping",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a first order conservative remapping "
  "of fields between grids in spherical coordinates.",
  "    The operators in this module uses code from the YAC software package to "
  "compute the conservative remapping weights.",
  "    For a detailed description of the interpolation method see YAC.",
  "    The interpolation method is completely general and can be used for any "
  "grid on a sphere.",
  "    The search algorithm for the conservative remapping requires that no "
  "grid cell occurs more than once. ",
  "",
  "OPERATORS",
  "    remapycon  First order conservative remapping",
  "               Performs a first order conservative remapping on all input "
  "fields.",
  "    genycon    Generate 1st order conservative remap weights",
  "               Generates first order conservative remapping weights for the "
  "first input field and",
  "               writes the result to a file. The format of this file is "
  "NetCDF following the SCRIP convention.",
  "               Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    CDO_REMAP_NORM",
  "        This variable is used to choose the normalization of the "
  "conservative interpolation. ",
  "        By default CDO_REMAP_NORM is set to 'fracarea'. 'fracarea' uses the "
  "sum of the",
  "        non-masked source cell intersected areas to normalize each target "
  "cell field value.",
  "        This results in a reasonable flux value but the flux is not locally "
  "conserved.",
  "        The option 'destarea' uses the total target cell area to normalize "
  "each target cell",
  "        field value. Local flux conservation is ensured, but unreasonable "
  "flux values may result.",
  "    REMAP_AREA_MIN",
  "        This variable is used to set the minimum destination area fraction. "
  "The default",
  "        of this variable is 0.0.",
};

std::vector<std::string> RemapconHelp = {
  "NAME",
  "    remapcon, gencon - First order conservative remapping",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a first order conservative remapping "
  "of fields between grids in spherical coordinates.",
  "    The interpolation is based on an adapted SCRIP library version. ",
  "    For a detailed description of the interpolation method see SCRIP.",
  "    The interpolation method is completely general and can be used for any "
  "grid on a sphere.",
  "    The search algorithm for the conservative remapping requires that no "
  "grid cell occurs more than once. ",
  "",
  "OPERATORS",
  "    remapcon  First order conservative remapping",
  "              Performs a first order conservative remapping on all input "
  "fields.",
  "    gencon    Generate 1st order conservative remap weights",
  "              Generates first order conservative remapping weights for the "
  "first input field and",
  "              writes the result to a file. The format of this file is "
  "NetCDF following the SCRIP convention.",
  "              Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    CDO_REMAP_NORM",
  "        This variable is used to choose the normalization of the "
  "conservative interpolation. ",
  "        By default CDO_REMAP_NORM is set to 'fracarea'. 'fracarea' uses the "
  "sum of the",
  "        non-masked source cell intersected areas to normalize each target "
  "cell field value.",
  "        This results in a reasonable flux value but the flux is not locally "
  "conserved.",
  "        The option 'destarea' uses the total target cell area to normalize "
  "each target cell",
  "        field value. Local flux conservation is ensured, but unreasonable "
  "flux values may result.",
  "    REMAP_AREA_MIN",
  "        This variable is used to set the minimum destination area fraction. "
  "The default",
  "        of this variable is 0.0.",
  "",
  "NOTE",
  "    The SCRIP conservative remapping method doesn't work correctly for some "
  "grid combinations.",
  "    Please use remapycon or genycon in case of problems. ",
};

std::vector<std::string> Remapcon2Help = {
  "NAME",
  "    remapcon2, gencon2 - Second order conservative remapping",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a second order conservative "
  "remapping of fields between grids in spherical coordinates.",
  "    The interpolation is based on an adapted SCRIP library version. ",
  "    For a detailed description of the interpolation method see SCRIP.",
  "    The interpolation method is completely general and can be used for any "
  "grid on a sphere.",
  "    The search algorithm for the conservative remapping requires that no "
  "grid cell occurs more than once. ",
  "",
  "OPERATORS",
  "    remapcon2  Second order conservative remapping",
  "               Performs a second order conservative remapping on all input "
  "fields.",
  "    gencon2    Generate 2nd order conservative remap weights",
  "               Generates second order conservative remapping weights for "
  "the first input field and",
  "               writes the result to a file. The format of this file is "
  "NetCDF following the SCRIP convention.",
  "               Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    CDO_REMAP_NORM",
  "        This variable is used to choose the normalization of the "
  "conservative interpolation. ",
  "        By default CDO_REMAP_NORM is set to 'fracarea'. 'fracarea' uses the "
  "sum of the",
  "        non-masked source cell intersected areas to normalize each target "
  "cell field value.",
  "        This results in a reasonable flux value but the flux is not locally "
  "conserved.",
  "        The option 'destarea' uses the total target cell area to normalize "
  "each target cell",
  "        field value. Local flux conservation is ensured, but unreasonable "
  "flux values may result.",
  "    REMAP_AREA_MIN",
  "        This variable is used to set the minimum destination area fraction. "
  "The default",
  "        of this variable is 0.0.",
  "",
  "NOTE",
  "    The SCRIP conservative remapping method doesn't work correctly for some "
  "grid combinations.",
};

std::vector<std::string> RemaplafHelp = {
  "NAME",
  "    remaplaf, genlaf - Largest area fraction remapping",
  "",
  "SYNOPSIS",
  "    <operator>,grid  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains operators for a largest area fraction remapping of "
  "fields between grids in spherical coordinates.",
  "    The operators in this module uses code from the YAC software package to "
  "compute the largest area fraction.",
  "    For a detailed description of the interpolation method see YAC.",
  "    The interpolation method is completely general and can be used for any "
  "grid on a sphere.",
  "    The search algorithm for this remapping method requires that no grid "
  "cell occurs more than once. ",
  "",
  "OPERATORS",
  "    remaplaf  Largest area fraction remapping",
  "              Performs a largest area fraction remapping on all input "
  "fields.",
  "    genlaf    Generate largest area fraction remap weights",
  "              Generates largest area fraction remapping weights for the "
  "first input field and",
  "              writes the result to a file. The format of this file is "
  "NetCDF following the SCRIP convention.",
  "              Use the operator remap to apply this remapping weights to a "
  "data file with the same source grid.",
  "",
  "PARAMETER",
  "    grid  STRING  Target grid description file or name",
  "",
  "ENVIRONMENT",
  "    REMAP_AREA_MIN",
  "        This variable is used to set the minimum destination area fraction. "
  "The default",
  "        of this variable is 0.0.",
};

std::vector<std::string> RemapHelp = {
  "NAME",
  "    remap - Grid remapping",
  "",
  "SYNOPSIS",
  "    remap,grid,weights  infile outfile",
  "",
  "DESCRIPTION",
  "    Interpolation between different horizontal grids can be a very "
  "time-consuming ",
  "    process. Especially if the data are on an unstructured and/or a large "
  "grid. ",
  "    In this case the interpolation process can be split into two parts.",
  "    Firstly the generation of the interpolation weights, which is the most "
  "time-consuming part.",
  "    These interpolation weights can be reused for every remapping process "
  "with the operator remap.",
  "    This operator remaps all input fields to a new horizontal grid. The "
  "remap type and ",
  "    the interpolation weights of one input grid are read from a NetCDF "
  "file. More weights ",
  "    are computed if the input fields are on different grids. The NetCDF "
  "file with the ",
  "    weights should follow the SCRIP convention. Normally these weights come "
  "from a previous",
  "    call to one of the genXXX operators (e.g. genbil) or were created by "
  "the original SCRIP package.",
  "",
  "PARAMETER",
  "    grid     STRING  Target grid description file or name",
  "    weights  STRING  Interpolation weights (SCRIP NetCDF file)",
  "",
  "ENVIRONMENT",
  "    CDO_REMAP_NORM       ",
  "        This variable is used to choose the normalization of the "
  "conservative interpolation. ",
  "        By default CDO_REMAP_NORM is set to 'fracarea'. 'fracarea' uses the "
  "sum of the",
  "        non-masked source cell intersected areas to normalize each target "
  "cell field value.",
  "        This results in a reasonable flux value but the flux is not locally "
  "conserved.",
  "        The option 'destarea' uses the total target cell area to normalize "
  "each target cell",
  "        field value. Local flux conservation is ensured, but unreasonable "
  "flux values may result.",
  "    REMAP_EXTRAPOLATE    ",
  "        This variable is used to switch the extrapolation feature 'on' or "
  "'off'.",
  "        By default the extrapolation is enabled for remapdis, remapnn and "
  "for circular grids.",
  "    REMAP_AREA_MIN       ",
  "        This variable is used to set the minimum destination area fraction. "
  "The default",
  "        of this variable is 0.0.",
  "    CDO_GRIDSEARCH_RADIUS",
  "        Grid search radius in degree, default 180 degree.",
};

std::vector<std::string> RemapetaHelp = {
  "NAME",
  "    remapeta - Remap vertical hybrid level",
  "",
  "SYNOPSIS",
  "    remapeta,vct[,oro]  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator interpolates between different vertical hybrid levels. "
  "This include the preparation ",
  "    of consistent data for the free atmosphere. The procedure for the "
  "vertical interpolation is based ",
  "    on the HIRLAM scheme and was adapted from INTERA.",
  "    The vertical interpolation is based on the vertical integration of the "
  "hydrostatic equation with ",
  "    few adjustments. The basic tasks are the following one:",
  "    - at first integration of hydrostatic equation",
  "    - extrapolation of surface pressure",
  "    - Planetary Boundary-Layer (PBL) proutfile interpolation",
  "    - interpolation in free atmosphere",
  "    - merging of both proutfiles",
  "    - final surface pressure correction",
  "    ",
  "    The vertical interpolation corrects the surface pressure. This is "
  "simply a cut-off or an addition ",
  "    of air mass. This mass correction should not influence the geostrophic "
  "velocity field in the middle ",
  "    troposhere. Therefore the total mass above a given reference level is "
  "conserved. As reference level",
  "    the geopotential height of the 400 hPa level is used. Near the surface "
  "the correction can affect ",
  "    the vertical structure of the PBL. Therefore the interpolation is done "
  "using the potential temperature. ",
  "    But in the free atmosphere above a certain n (n=0.8 defining the top of "
  "the PBL) the interpolation ",
  "    is done linearly. After the interpolation both proutfiles are merged. "
  "With the resulting ",
  "    temperature/pressure correction the hydrostatic equation is integrated "
  "again and adjusted to the ",
  "    reference level finding the final surface pressure correction. A more "
  "detailed description of",
  "    the interpolation can be found in INTERA. This operator requires all "
  "variables on the same horizontal grid.",
  "",
  "PARAMETER",
  "    vct  STRING  File name of an ASCII dataset with the vertical coordinate "
  "table",
  "    oro  STRING  File name with the orography (surf. geopotential) of the "
  "target dataset (optional)",
  "",
  "ENVIRONMENT",
  "    REMAPETA_PTOP",
  "        Sets the minimum pressure level for condensation.",
  "        Above this level the humidity is set to the constant 1.E-6.",
  "        The default value is 0 Pa.",
  "",
  "NOTE",
  "    The code numbers or the variable names of the required parameter have "
  "to follow the ECHAM convention.",
  "    Presently, the vertical coordinate definition of a NetCDF file has also "
  "to follow the ECHAM convention.",
  "    This means:",
  "    - the dimension of the full level coordinate and the corresponding "
  "variable is called mlev,",
  "    - the dimension of the half level coordinate and the corresponding "
  "variable is called ilev (ilev must have one element more than mlev)",
  "    - the hybrid vertical coefficient a is given in units of Pa and called "
  "hyai (hyam for level midpoints)",
  "    - the hybrid vertical coefficient b is given in units of 1 and called "
  "hybi (hybm for level midpoints)",
  "    - the mlev variable has a borders attribute containing the character "
  "string 'ilev'",
  "    ",
  "    Use the sinfo command to test if your vertical coordinate system is "
  "recognized as hybrid system.",
  "    ",
  "    In case remapeta complains about not finding any data on hybrid model "
  "levels you may wish",
  "    to use the setzaxis command to generate a zaxis description which "
  "conforms to the ECHAM convention.",
  "    See section \"1.4 Z-axis description\" for an example how to define a "
  "hybrid Z-axis.",
};

std::vector<std::string> VertintmlHelp = {
  "NAME",
  "    ml2pl, ml2hl - Vertical interpolation",
  "",
  "SYNOPSIS",
  "    ml2pl,plevels  infile outfile",
  "    ml2hl,hlevels  infile outfile",
  "",
  "DESCRIPTION",
  "    Interpolate 3D variables on hybrid sigma pressure level to pressure or "
  "height levels.",
  "    The input file should contain the log. surface pressure or the surface "
  "pressure.",
  "    To extrapolate the temperature, the surface geopotential is also "
  "needed.",
  "    The pressure, temperature, and surface geopotential are identified by "
  "their GRIB1 code number",
  "    or NetCDF CF standard name.",
  "    Supported parameter tables are: WMO standard table number 2 and ECMWF "
  "local table number 128.",
  "    Use the alias  ml2plx/ml2hlx or the environment variable EXTRAPOLATE",
  "    to extrapolate missing values. This operator requires all variables on "
  "the same horizontal grid.",
  "    ",
  "",
  "OPERATORS",
  "    ml2pl  Model to pressure level interpolation",
  "           Interpolates 3D variables on hybrid sigma pressure level to "
  "pressure level.",
  "    ml2hl  Model to height level interpolation",
  "           Interpolates 3D variables on hybrid sigma pressure level to "
  "height level.",
  "           The procedure is the same as for the operator ml2pl except for",
  "           the pressure levels being calculated from the heights by:",
  "           plevel = 101325*exp(hlevel/-7000)",
  "",
  "PARAMETER",
  "    plevels  FLOAT  Pressure levels in pascal",
  "    hlevels  FLOAT  Height levels in meter (max level: 65535 m)",
  "",
  "ENVIRONMENT",
  "    EXTRAPOLATE",
  "        If set to 1 extrapolate missing values.",
};

std::vector<std::string> VertintapHelp = {
  "NAME",
  "    ap2pl, ap2hl - Vertical interpolation",
  "",
  "SYNOPSIS",
  "    ap2pl,plevels  infile outfile",
  "    ap2hl,hlevels  infile outfile",
  "",
  "DESCRIPTION",
  "    Interpolate 3D variables on hybrid sigma height coordinates to pressure "
  "or height levels.",
  "    The input file must contain the 3D air pressure. The air pressure is "
  "identified",
  "    by the NetCDF CF standard name air_pressure.",
  "    Use the alias  ap2plx/ap2hlx or the environment variable EXTRAPOLATE",
  "    to extrapolate missing values. This operator requires all variables on "
  "the same horizontal grid.",
  "",
  "OPERATORS",
  "    ap2pl  Air pressure to pressure level interpolation",
  "           Interpolates 3D variables on hybrid sigma height coordinates to "
  "pressure level.",
  "    ap2hl  Air pressure to height level interpolation",
  "           Interpolates 3D variables on hybrid sigma height coordinates to "
  "height level.",
  "           The procedure is the same as for the operator ap2pl except for",
  "           the pressure levels being calculated from the heights by:",
  "           plevel = 101325*exp(hlevel/-7000)",
  "",
  "PARAMETER",
  "    plevels  FLOAT  Pressure levels in pascal",
  "    hlevels  FLOAT  Height levels in meter (max level: 65535 m)",
  "",
  "ENVIRONMENT",
  "    EXTRAPOLATE",
  "        If set to 1 extrapolate missing values.",
  "",
  "NOTE",
  "    This is a specific implementation for NetCDF files from the ICON model, "
  "it may not work with data from other sources.",
};

std::vector<std::string> IntlevelHelp = {
  "NAME",
  "    intlevel - Linear level interpolation",
  "",
  "SYNOPSIS",
  "    intlevel,levels  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator performs a linear vertical interpolation of non hybrid 3D "
  "variables.",
  "",
  "PARAMETER",
  "    levels  FLOAT  Target levels",
};

std::vector<std::string> Intlevel3dHelp = {
  "NAME",
  "    intlevel3d, intlevelx3d - ",
  "    Linear level interpolation from/to 3d vertical coordinates",
  "",
  "SYNOPSIS",
  "    <operator>,icoordinate  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    This operator performs a linear vertical interpolation of 3D variables "
  "fields",
  "    with given 3D vertical coordinates.",
  "",
  "OPERATORS",
  "    intlevel3d   Linear level interpolation onto a 3d vertical coordinate",
  "    intlevelx3d  like intlevel3d but with extrapolation",
  "",
  "PARAMETER",
  "    icoordinate  STRING  filename for vertical source coordinates variable",
  "    infile2      STRING  target vertical coordinate field (intlevel3d only)",
};

std::vector<std::string> InttimeHelp = {
  "NAME",
  "    inttime, intntime - Time interpolation",
  "",
  "SYNOPSIS",
  "    inttime,date,time[,inc]  infile outfile",
  "    intntime,n  infile outfile",
  "",
  "DESCRIPTION",
  "    This module performs linear interpolation between timesteps.",
  "",
  "OPERATORS",
  "    inttime   Interpolation between timesteps",
  "              This operator creates a new dataset by linear interpolation "
  "between timesteps.",
  "              The user has to define the start date/time with an optional "
  "increment.",
  "    intntime  Interpolation between timesteps",
  "              This operator performs linear interpolation between "
  "timesteps.",
  "              The user has to define the number of timesteps from one "
  "timestep to the next.",
  "",
  "PARAMETER",
  "    date  STRING  Start date (format YYYY-MM-DD)",
  "    time  STRING  Start time (format hh:mm:ss)",
  "    inc   STRING  Optional increment (seconds, minutes, hours, days, "
  "months, years) [default: 0hour]",
  "    n     INTEGER Number of timesteps from one timestep to the next",
};

std::vector<std::string> IntyearHelp = {
  "NAME",
  "    intyear - Year interpolation",
  "",
  "SYNOPSIS",
  "    intyear,years  infile1 infile2 obase",
  "",
  "DESCRIPTION",
  "    This operator performs linear interpolation between two years, timestep "
  "by timestep.",
  "    The input files need to have the same structure with the same "
  "variables.",
  "    The output files will be named <obase><yyyy><suffix> where yyyy will be "
  "the year and ",
  "    suffix is the filename extension derived from the file format.",
  "",
  "PARAMETER",
  "    years  INTEGER  Comma separated list of years",
  "",
  "ENVIRONMENT",
  "    CDO_FILE_SUFFIX",
  "        Set the default file suffix. This suffix will be added to the "
  "output file ",
  "        names instead of the filename extension derived from the file "
  "format. ",
  "        Set this variable to NULL to disable the adding of a file suffix.",
  "",
  "NOTE",
  "    This operator needs to open all output files simultaneously.",
  "    The maximum number of open files depends on the operating system!",
};

std::vector<std::string> SpectralHelp = {
  "NAME",
  "    sp2gp, sp2gpl, gp2sp, gp2spl, sp2sp - Spectral transformation",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "    sp2sp,trunc  infile outfile",
  "",
  "DESCRIPTION",
  "    This module transforms fields on a global regular Gaussian grids to "
  "spectral coefficients and vice versa.",
  "    Missing values are not supported.",
  "",
  "OPERATORS",
  "    sp2gp   Spectral to gridpoint",
  "            Convert all fields with spectral coefficients to a global "
  "regular Gaussian grid. The number of ",
  "            latitudes of the resulting Gaussian grid is calculated from the "
  "triangular truncation by:",
  "            ",
  "               nlat = NINT((trunc*3 + 1.)/2.)",
  "    sp2gpl  Spectral to gridpoint (linear)",
  "            Convert all fields with spectral coefficients to a global "
  "regular Gaussian grid. The number of ",
  "            latitudes of the resulting Gaussian grid is calculated from the "
  "triangular truncation by:",
  "            ",
  "               nlat = NINT((trunc*2 + 1.)/2.)",
  "            ",
  "            Use this operator to convert ERA40 data e.g. from TL159 to N80.",
  "    gp2sp   Gridpoint to spectral",
  "            Convert all Gaussian gridpoint fields to spectral coefficients. "
  "The triangular truncation ",
  "            of the resulting spherical harmonics is calculated from the "
  "number of latitudes by:",
  "            ",
  "               trunc = (nlat*2 - 1) / 3",
  "    gp2spl  Gridpoint to spectral (linear)",
  "            Convert all Gaussian gridpoint fields to spectral coefficients. "
  "The triangular truncation ",
  "            of the resulting spherical harmonics is calculated from the "
  "number of latitudes by:",
  "            ",
  "               trunc = (nlat*2 - 1) / 2",
  "            ",
  "            Use this operator to convert ERA40 data e.g. from N80 to TL159 "
  "instead of T106.",
  "    sp2sp   Spectral to spectral",
  "            Change the triangular truncation of all spectral fields. The "
  "operator performs downward ",
  "            conversion by cutting the resolution. Upward conversions are "
  "achieved by filling in zeros.",
  "",
  "PARAMETER",
  "    trunc  INTEGER  New spectral resolution",
};

std::vector<std::string> WindHelp = {
  "NAME",
  "    dv2uv, dv2uvl, uv2dv, uv2dvl, dv2ps - Wind transformation",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module converts relative divergence and vorticity to U and V wind "
  "and vice versa.",
  "    Divergence and vorticity are spherical harmonic coefficients in "
  "spectral space and",
  "    U and V are on a global regular Gaussian grid. The Gaussian latitudes "
  "need to be ordered from",
  "    north to south. Missing values are not supported.",
  "",
  "OPERATORS",
  "    dv2uv   Divergence and vorticity to U and V wind",
  "            Calculate U and V wind on a Gaussian grid from spherical "
  "harmonic ",
  "            coefficients of relative divergence and vorticity. The "
  "divergence and vorticity ",
  "            need to have the names sd and svo or code numbers 155 and 138.",
  "            The number of latitudes of the resulting Gaussian grid is "
  "calculated ",
  "            from the triangular truncation by:",
  "            ",
  "               nlat = NINT((trunc*3 + 1.)/2.)",
  "    dv2uvl  Divergence and vorticity to U and V wind (linear)",
  "            Calculate U and V wind on a Gaussian grid from spherical "
  "harmonic ",
  "            coefficients of relative divergence and vorticity. The "
  "divergence and vorticity ",
  "            need to have the names sd and svo or code numbers 155 and 138.",
  "            The number of latitudes of the resulting Gaussian grid is "
  "calculated ",
  "            from the triangular truncation by:",
  "            ",
  "               nlat = NINT((trunc*2 + 1.)/2.)",
  "    uv2dv   U and V wind to divergence and vorticity",
  "            Calculate spherical harmonic coefficients of relative "
  "divergence and vorticity",
  "            from U and V wind. The U and V wind need to be on a Gaussian "
  "grid and need to have the ",
  "            names u and v or the code numbers 131 and 132.",
  "            The triangular truncation of the resulting spherical harmonics",
  "            is calculated from the number of latitudes by:",
  "            ",
  "               trunc = (nlat*2 - 1) / 3",
  "    uv2dvl  U and V wind to divergence and vorticity (linear)",
  "            Calculate spherical harmonic coefficients of relative "
  "divergence and vorticity",
  "            from U and V wind. The U and V wind need to be on a Gaussian "
  "grid and need to have the ",
  "            names u and v or the code numbers 131 and 132.",
  "            The triangular truncation of the resulting spherical harmonics",
  "            is calculated from the number of latitudes by:",
  "            ",
  "               trunc = (nlat*2 - 1) / 2",
  "    dv2ps   D and V to velocity potential and stream function",
  "            Calculate spherical harmonic coefficients of velocity potential "
  "and stream function from ",
  "            spherical harmonic coefficients of relative divergence and "
  "vorticity. The divergence and ",
  "            vorticity need to have the names sd and svo or code numbers 155 "
  "and 138.",
};

std::vector<std::string> ImportbinaryHelp = {
  "NAME",
  "    import_binary - Import binary data sets",
  "",
  "SYNOPSIS",
  "    import_binary  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator imports gridded binary data sets via a GrADS data "
  "descriptor file.",
  "    The GrADS data descriptor file contains a complete description of the "
  "binary data as well ",
  "    as instructions on where to find the data and how to read it. The "
  "descriptor file is an ASCII ",
  "    file that can be created easily with a text editor. The general "
  "contents of a gridded data ",
  "    descriptor file are as follows:",
  "    - Filename for the binary data",
  "    - Missing or undefined data value",
  "    - Mapping between grid coordinates and world coordinates",
  "    - Description of variables in the binary data set ",
  "    ",
  "    A detailed description of the components of a GrADS data descriptor "
  "file can be found in GrADS.",
  "    Here is a list of the supported components:",
  "    BYTESWAPPED, CHSUB, DSET, ENDVARS, FILEHEADER, HEADERBYTES, OPTIONS, "
  "TDEF, TITLE, ",
  "    TRAILERBYTES, UNDEF, VARS, XDEF, XYHEADER, YDEF, ZDEF",
  "",
  "NOTE",
  "    Only 32-bit IEEE floats are supported for standard binary files!",
};

std::vector<std::string> ImportcmsafHelp = {
  "NAME",
  "    import_cmsaf - Import CM-SAF HDF5 files",
  "",
  "SYNOPSIS",
  "    import_cmsaf  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator imports gridded CM-SAF (Satellite Application Facility on "
  "Climate Monitoring)",
  "    HDF5 files. CM-SAF exploits data from polar-orbiting and geostationary "
  "satellites in order ",
  "    to provide climate monitoring products of the following parameters: ",
  "    ",
  "    Cloud parameters: cloud fraction (CFC), cloud type (CTY), cloud phase "
  "(CPH), ",
  "                      cloud top height, pressure and temperature "
  "(CTH,CTP,CTT), ",
  "                      cloud optical thickness (COT), cloud water path "
  "(CWP).",
  "    ",
  "    Surface radiation components: Surface albedo (SAL); surface incoming "
  "(SIS) ",
  "                      and net (SNS) shortwave radiation; surface downward "
  "(SDL) ",
  "                      and outgoing (SOL) longwave radiation, surface net "
  "longwave ",
  "                      radiation (SNL) and surface radiation budget (SRB).",
  "    ",
  "    Top-of-atmosphere radiation components: Incoming (TIS) and reflected "
  "(TRS) ",
  "                      solar radiative flux at top-of-atmosphere. Emitted "
  "thermal ",
  "                      radiative flux at top-of-atmosphere (TET).",
  "    ",
  "    Water vapour:     Vertically integrated water vapour (HTW), layered "
  "vertically ",
  "                      integrated water vapour and layer mean temperature "
  "and relative ",
  "                      humidity for 5 layers (HLW), temperature and mixing "
  "ratio at ",
  "                      6 pressure levels. ",
  "    ",
  "    Daily and monthly mean products can be ordered via the CM-SAF web page "
  "(www.cmsaf.eu). ",
  "    Products with higher spatial and temporal resolution, i.e. "
  "instantaneous swath-based products,",
  "    are available on request (contact.cmsaf@dwd.de). All products are "
  "distributed free-of-charge.",
  "    More information on the data is available on the CM-SAF homepage "
  "(www.cmsaf.eu).",
  "    ",
  "    Daily and monthly mean products are provided in equal-area projections. "
  "CDO reads the ",
  "    projection parameters from the metadata in the HDF5-headers in order to "
  "allow spatial ",
  "    operations like remapping. For spatial operations with instantaneous "
  "products on original ",
  "    satellite projection, additional files with arrays of latitudes and "
  "longitudes are needed.",
  "    These can be obtained from CM-SAF together with the data.",
  "    ",
  "",
  "NOTE",
  "    To use this operator, it is necessary to build CDO with HDF5 support "
  "(version 1.6 or higher).",
  "    The PROJ.4 library (version 4.6 or higher) is needed for full support "
  "of the remapping",
  "    functionality. ",
};

std::vector<std::string> ImportamsrHelp = {
  "NAME",
  "    import_amsr - Import AMSR binary files",
  "",
  "SYNOPSIS",
  "    import_amsr  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator imports gridded binary AMSR (Advanced Microwave Scanning "
  "Radiometer) data.",
  "    The binary data files are available from the AMSR ftp site "
  "(ftp://ftp.ssmi.com/amsre).",
  "    Each file consists of twelve (daily) or five (averaged) 0.25 x 0.25 "
  "degree ",
  "    grid (1440,720) byte maps. For daily files, six daytime maps in the "
  "following",
  "    order, Time (UTC), Sea Surface Temperature (SST), 10 meter Surface Wind "
  "Speed (WSPD),",
  "    Atmospheric Water Vapor (VAPOR), Cloud Liquid Water (CLOUD), and Rain "
  "Rate (RAIN), ",
  "    are followed by six nighttime maps in the same order. Time-Averaged "
  "files contain ",
  "    just the geophysical layers in the same order [SST, WSPD, VAPOR, CLOUD, "
  "RAIN].",
  "    More information to the data is available on the AMSR homepage "
  "http://www.remss.com/amsr.",
};

std::vector<std::string> InputHelp = {
  "NAME",
  "    input, inputsrv, inputext - Formatted input",
  "",
  "SYNOPSIS",
  "    input,grid[,zaxis]  outfile",
  "    inputsrv  outfile",
  "    inputext  outfile",
  "",
  "DESCRIPTION",
  "    This module reads time series of one 2D variable from standard input.",
  "    All input fields need to have the same horizontal grid. The format of "
  "the ",
  "    input depends on the chosen operator.",
  "",
  "OPERATORS",
  "    input     ASCII input",
  "              Reads fields with ASCII numbers from standard input and "
  "stores them",
  "              in outfile. The numbers read are exactly that ones which are "
  "written ",
  "              out by the output operator.",
  "    inputsrv  SERVICE ASCII input",
  "              Reads fields with ASCII numbers from standard input and "
  "stores them ",
  "              in outfile. Each field should have a header of 8 integers "
  "(SERVICE likely).",
  "              The numbers that are read are exactly that ones which are "
  "written out by ",
  "              the outputsrv operator.",
  "    inputext  EXTRA ASCII input",
  "              Read fields with ASCII numbers from standard input and stores "
  "them ",
  "              in outfile. Each field should have header of 4 integers "
  "(EXTRA likely).",
  "              The numbers read are exactly that ones which are written out "
  "by ",
  "              the outputext operator.",
  "",
  "PARAMETER",
  "    grid   STRING  Grid description file or name",
  "    zaxis  STRING  Z-axis description file",
};

std::vector<std::string> OutputHelp = {
  "NAME",
  "    output, outputf, outputint, outputsrv, outputext - Formatted output",
  "",
  "SYNOPSIS",
  "    output  infiles",
  "    outputf,format[,nelem]  infiles",
  "    outputint  infiles",
  "    outputsrv  infiles",
  "    outputext  infiles",
  "",
  "DESCRIPTION",
  "    This module prints all values of all input datasets to standard output.",
  "    All input fields need to have the same horizontal grid. All input "
  "files ",
  "    need to have the same structure with the same variables.",
  "    The format of the output depends on the chosen operator.",
  "",
  "OPERATORS",
  "    output     ASCII output",
  "               Prints all values to standard output.",
  "               Each row has 6 elements with the C-style format \"%13.6g\".",
  "    outputf    Formatted output",
  "               Prints all values to standard output.",
  "               The format and number of elements for each row have to be "
  "specified by the parameters",
  "               format and nelem. The default for nelem is 1.",
  "    outputint  Integer output",
  "               Prints all values rounded to the nearest integer to standard "
  "output.",
  "    outputsrv  SERVICE ASCII output",
  "               Prints all values to standard output.",
  "               Each field with a header of 8 integers (SERVICE likely).",
  "    outputext  EXTRA ASCII output",
  "               Prints all values to standard output.",
  "               Each field with a header of 4 integers (EXTRA likely).",
  "",
  "PARAMETER",
  "    format  STRING  C-style format for one element (e.g. %13.6g)",
  "    nelem   INTEGER Number of elements for each row (default: nelem = 1)",
};

std::vector<std::string> OutputtabHelp = {
  "NAME",
  "    outputtab - Table output",
  "",
  "SYNOPSIS",
  "    outputtab,params  infiles outfile",
  "",
  "DESCRIPTION",
  "    This operator prints a table of all input datasets to standard output.",
  "    infiles is an arbitrary number of input files. All input files need to "
  "have ",
  "    the same structure with the same variables on different timesteps.",
  "    All input fields need to have the same horizontal grid.",
  "    ",
  "    The contents of the table depends on the chosen paramters. The format "
  "of each table",
  "    parameter is keyname[:len]. len is the optional length of a table "
  "entry.  ",
  "    Here is a list of all valid keynames:",
  "    ",
  "     Keyname    & Type    & Description      ",
  "     value      & FLOAT   & Value of the variable [len:8]",
  "     name       & STRING  & Name of the variable [len:8]",
  "     param      & STRING  & Parameter ID (GRIB1: code[.tabnum]; GRIB2: "
  "num[.cat[.dis]]) [len:11]",
  "     code       & INTEGER & Code number [len:4]",
  "     lon        & FLOAT   & Longitude coordinate [len:6]",
  "     lat        & FLOAT   & Latitude coordinate [len:6]",
  "     lev        & FLOAT   & Vertical level [len:6]",
  "     xind       & INTEGER & Grid x index [len:4]",
  "     yind       & INTEGER & Grid y index [len:4]",
  "     timestep   & INTEGER & Timestep number [len:6]",
  "     date       & STRING  & Date (format YYYY-MM-DD) [len:10]",
  "     time       & STRING  & Time (format hh:mm:ss) [len:8]",
  "     year       & INTEGER & Year [len:5]",
  "     month      & INTEGER & Month [len:2]",
  "     day        & INTEGER & Day [len:2]",
  "     nohead     & INTEGER & Disable output of header line",
  "",
  "PARAMETER",
  "    params  STRING   Comma separated list of keynames, one for each column "
  "of the table",
};

std::vector<std::string> OutputgmtHelp = {
  "NAME",
  "    gmtxyz, gmtcells - GMT output",
  "",
  "SYNOPSIS",
  "    <operator>  infile",
  "",
  "DESCRIPTION",
  "    This module prints the first field of the input dataset to standard "
  "output.",
  "    The output can be used to generate 2D Lon/Lat plots with GMT.",
  "    The format of the output depends on the chosen operator.",
  "",
  "OPERATORS",
  "    gmtxyz    GMT xyz format",
  "              The operator exports the first field to the GMT xyz ASCII "
  "format.",
  "              The output can be used to create contour plots with the GMT "
  "module pscontour.",
  "    gmtcells  GMT multiple segment format",
  "              The operator exports the first field to the GMT multiple "
  "segment ASCII format.",
  "              The output can be used to create shaded gridfill plots with "
  "the GMT module psxy.",
};

std::vector<std::string> GradsdesHelp = {
  "NAME",
  "    gradsdes - GrADS data descriptor file",
  "",
  "SYNOPSIS",
  "    gradsdes[,mapversion]  infile",
  "",
  "DESCRIPTION",
  "    Creates a GrADS data descriptor file. Supported file formats are GRIB1, "
  "NetCDF, SERVICE, ",
  "    EXTRA and IEG. For GRIB1 files the GrADS map file is also generated. "
  "For SERVICE and EXTRA",
  "    files the grid have to be specified with the CDO option '-g <grid>'. "
  "This module takes infile",
  "    in order to create filenames for the descriptor (infile.ctl) and the "
  "map (infile.gmp) file.",
  "",
  "PARAMETER",
  "    mapversion  INTEGER  Format version of the GrADS map file for GRIB1 "
  "datasets. Use 1 for a machine",
  "                specific version 1 GrADS map file, 2 for a machine "
  "independent version 2 GrADS map file",
  "                and 4 to support GRIB files >2GB. ",
  "                A version 2 map file can be used only with GrADS version "
  "1.8 or newer.",
  "                A version 4 map file can be used only with GrADS version "
  "2.0 or newer.",
  "                The default is 4 for files >2GB, otherwise 2.",
};

std::vector<std::string> AfterburnerHelp = {
  "NAME",
  "    after - ECHAM standard post processor",
  "",
  "SYNOPSIS",
  "    after[,vct]  infiles outfile",
  "",
  "DESCRIPTION",
  "    The \"afterburner\" is the standard post processor for ECHAM data which "
  "provides the following operations:",
  "    ",
  "    - Extract specified variables and levels",
  "    - Compute derived variables",
  "    - Transform spectral data to Gaussian grid representation",
  "    - Vertical interpolation to pressure levels",
  "    - Compute temporal means",
  "    ",
  "    This operator reads selection parameters as namelist from stdin.",
  "    Use the UNIX redirection \"<namelistfile\" to read the namelist from "
  "file.",
  "",
  "NAMELIST",
  "    Namelist parameter and there defaults:",
  "      TYPE=0, CODE=-1, LEVEL=-1, INTERVAL=0, MEAN=0, EXTRAPOLATE=0",
  "    ",
  "    TYPE controls the transformation and vertical interpolation. "
  "Transforming spectral data to Gaussian grid",
  "    representation and vertical interpolation to pressure levels are "
  "performed in a chain of steps.",
  "    The TYPE parameter may be used to stop the chain at a certain step. "
  "Valid values are:",
  "    ",
  "      TYPE  =  0 : Hybrid   level spectral coefficients",
  "      TYPE  = 10 : Hybrid   level fourier  coefficients",
  "      TYPE  = 11 : Hybrid   level zonal mean sections",
  "      TYPE  = 20 : Hybrid   level gauss grids",
  "      TYPE  = 30 : Pressure level gauss grids",
  "      TYPE  = 40 : Pressure level fourier  coefficients",
  "      TYPE  = 41 : Pressure level zonal mean sections",
  "      TYPE  = 50 : Pressure level spectral coefficients",
  "      TYPE  = 60 : Pressure level fourier  coefficients",
  "      TYPE  = 61 : Pressure level zonal mean sections",
  "      TYPE  = 70 : Pressure level gauss grids",
  "    ",
  "    Vorticity, divergence, streamfunction and velocity potential need "
  "special treatment in the vertical transformation.",
  "    They are not available as types 30, 40 and 41. If you select one of "
  "these combinations, type is automatically",
  "    switched to the equivalent types 70, 60 and 61. The type of all other "
  "variables will be switched too, because ",
  "    the type is a global parameter.",
  "    ",
  "    CODE selects the variables by the ECHAM GRIB1 code number (1-255). The "
  "default value -1 processes all detected codes.",
  "    Derived variables computed by the afterburner:",
  "    ",
  "    Code  & Name      & Longname                       & Units & Level      "
  " & Needed Codes",
  "     34   & low_cld   & low cloud                      &       & single     "
  " & 223 on modellevel  ",
  "     35   & mid_cld   & mid cloud                      &       & single     "
  " & 223 on modellevel  ",
  "     36   & hih_cld   & high cloud                     &       & single     "
  " & 223 on modellevel  ",
  "     131  & u         & u-velocity                     & m/s   & atm "
  "(ml+pl) & 138, 155           ",
  "     132  & v         & v-velocity                     & m/s   & atm "
  "(ml+pl) & 138, 155           ",
  "     135  & omega     & vertical velocity              & Pa/s  & atm "
  "(ml+pl) & 138, 152, 155      ",
  "     148  & stream    & streamfunction                 & m^2/s & atm "
  "(ml+pl) & 131, 132           ",
  "     149  & velopot   & velocity potential             & m^2/s & atm "
  "(ml+pl) & 131, 132           ",
  "     151  & slp       & mean sea level pressure        & Pa    & surface    "
  " & 129, 130, 152       ",
  "     156  & geopoth   & geopotential height            & m     & atm "
  "(ml+pl) & 129, 130, 133, 152 ",
  "     157  & rhumidity & relative humidity              &       & atm "
  "(ml+pl) & 130, 133, 152      ",
  "     189  & sclfs     & surface solar cloud forcing    &       & surface    "
  " & 176-185            ",
  "     190  & tclfs     & surface thermal cloud forcing  &       & surface    "
  " & 177-186            ",
  "     191  & sclf0     & top solar cloud forcing        &       & surface    "
  " & 178-187             ",
  "     192  & tclf0     & top thermal cloud forcing      &       & surface    "
  " & 179-188            ",
  "     259  & windspeed & windspeed                      & m/s   & atm "
  "(ml+pl) & sqrt(u*u+v*v)      ",
  "     260  & precip    & total precipitation            &       & surface    "
  " & 142+143            ",
  "    ",
  "    LEVEL selects the hybrid or pressure levels. The allowed values depends "
  "on the parameter TYPE.",
  "    The default value -1 processes all detected levels.",
  "    ",
  "    INTERVAL selects the processing interval. The default value 0 process "
  "data on monthly intervals.",
  "    INTERVAL=1 sets the interval to daily.",
  "    ",
  "    MEAN=1 compute and write monthly or daily mean fields. The default "
  "value 0 writes out all timesteps.",
  "    ",
  "    EXTRAPOLATE=0 switch of the extrapolation of missing values during the "
  "interpolation from model to pressure",
  "    level (only available with MEAN=0 and TYPE=30). The default value 1 "
  "extrapolate missing values.",
  "    ",
  "    Possible combinations of TYPE, CODE and MEAN:",
  "    ",
  "          TYPE   & CODE                    & MEAN",
  "        0/10/11  & 130  temperature        &  0",
  "        0/10/11  & 131  u-velocity         &  0",
  "        0/10/11  & 132  v-velocity         &  0",
  "        0/10/11  & 133  specific humidity  &  0",
  "        0/10/11  & 138  vorticity          &  0",
  "        0/10/11  & 148  streamfunction     &  0",
  "        0/10/11  & 149  velocity potential &  0",
  "        0/10/11  & 152  LnPs               &  0",
  "        0/10/11  & 155  divergence         &  0",
  "         >11     & all codes               &  0/1",
  "",
  "PARAMETER",
  "    vct  STRING  File with VCT in ASCII format",
};

std::vector<std::string> FilterHelp = {
  "NAME",
  "    bandpass, lowpass, highpass - Time series filtering",
  "",
  "SYNOPSIS",
  "    bandpass,fmin,fmax  infile outfile",
  "    lowpass,fmax  infile outfile",
  "    highpass,fmin  infile outfile",
  "",
  "DESCRIPTION",
  "    This module takes the time series for each gridpoint in infile and "
  "(fast fourier) transforms it ",
  "    into the frequency domain. According to the particular operator and its "
  "parameters certain frequencies ",
  "    are filtered (set to zero) in the frequency domain and the spectrum is "
  "(inverse fast fourier) transformed ",
  "    back into the time domain.",
  "    To determine the frequency the time-axis of infile is used. (Data "
  "should have a constant time increment ",
  "    since this assumption applies for transformation. However, the time "
  "increment has to be different from zero.)",
  "    All frequencies given as parameter are interpreted per year. This is "
  "done by the assumption of a 365-day calendar. ",
  "    Consequently if you want to perform multiyear-filtering accurately you "
  "have to delete the 29th of February. ",
  "    If your infile has a 360 year calendar the frequency parameters fmin "
  "respectively fmax should be ",
  "    multiplied with a factor of 360/365 in order to obtain accurate "
  "results.  ",
  "    For the set up of a frequency filter the frequency parameters have to "
  "be adjusted to a frequency in the data. ",
  "    Here fmin is rounded down and fmax is always rounded up. Consequently "
  "it is possible to use bandpass with ",
  "    fmin=fmax without getting a zero-field for outfile. ",
  "    Hints for efficient usage: ",
  "    - to get reliable results the time-series has to be detrended (cdo "
  "detrend)",
  "    - the lowest frequency greater zero that can be contained in infile is "
  "1/(N*dT), ",
  "    - the greatest frequency is 1/(2dT) (Nyquist frequency),",
  "    with N the number of timesteps and dT the time increment of infile in "
  "years. ",
  "",
  "OPERATORS",
  "    bandpass  Bandpass filtering",
  "              Bandpass filtering (pass for frequencies between fmin and "
  "fmax).",
  "              Suppresses all variability outside the frequency range "
  "specified by [fmin,fmax].",
  "    lowpass   Lowpass filtering",
  "              Lowpass filtering (pass for frequencies lower than fmax).",
  "              Suppresses all variability with frequencies greater than "
  "fmax. ",
  "    highpass  Highpass filtering",
  "              Highpass filtering (pass for frequencies greater than fmin). ",
  "              Suppresses all variabilty with frequencies lower than fmin. ",
  "",
  "PARAMETER",
  "    fmin  FLOAT	Minimum frequency per year that passes the filter.",
  "    fmax  FLOAT	Maximum frequency per year that passes the filter.  ",
  "",
  "NOTE",
  "    For better performace of these operators use the CDO configure option "
  "--with-fftw3.",
};

std::vector<std::string> GridcellHelp = {
  "NAME",
  "    gridarea, gridweights - Grid cell quantities",
  "",
  "SYNOPSIS",
  "    <operator>  infile outfile",
  "",
  "DESCRIPTION",
  "    This module reads the grid cell area of the first grid from the input "
  "stream.",
  "    If the grid cell area is missing it will be computed from the ",
  "    grid description. Depending on the chosen operator the grid cell area "
  "or weights",
  "    are written to the output stream.",
  "",
  "OPERATORS",
  "    gridarea     Grid cell area",
  "                 Writes the grid cell area to the output stream. If the "
  "grid cell area have to",
  "                 be computed it is scaled with the earth radius to square "
  "meters.",
  "    gridweights  Grid cell weights",
  "                 Writes the grid cell area weights to the output stream.",
  "",
  "ENVIRONMENT",
  "    PLANET_RADIUS",
  "        This variable is used to scale the computed grid cell areas to "
  "square meters. ",
  "        By default PLANET_RADIUS is set to an earth radius of 6371000 "
  "meter.",
};

std::vector<std::string> SmoothHelp = {
  "NAME",
  "    smooth, smooth9 - Smooth grid points",
  "",
  "SYNOPSIS",
  "    smooth[,options]  infile outfile",
  "    smooth9  infile outfile",
  "",
  "DESCRIPTION",
  "    Smooth all grid points of a horizontal grid.",
  "    Options is a comma separated list of \"key=value\" pairs with optional "
  "parameters.",
  "",
  "OPERATORS",
  "    smooth   Smooth grid points",
  "             Performs a N point smoothing on all input fields. The number "
  "of points used depend",
  "             on the search radius (radius) and the maximum number of points "
  "(maxpoints).",
  "             Per default all points within the search radius of 1degree are "
  "used.",
  "             The weights for the points depend on the form of the curve and "
  "the distance.",
  "             The implemented form of the curve is linear with constant "
  "default weights of 0.25",
  "             at distance 0 (weight0) and at the search radius (weightR).",
  "    smooth9  9 point smoothing",
  "             Performs a 9 point smoothing on all fields with a "
  "quadrilateral curvilinear grid.",
  "             The result at each grid point is a weighted average of the "
  "grid point plus",
  "             the 8 surrounding points. The center point receives a weight "
  "of 1.0, the ",
  "             points at each side and above and below receive a weight of "
  "0.5, and corner ",
  "             points receive a weight of 0.3.",
  "             All 9 points are multiplied by their weights and summed, then "
  "divided by ",
  "             the total weight to obtain the smoothed value. Any missing "
  "data points are ",
  "             not included in the sum; points beyond the grid boundary are "
  "considered to ",
  "             be missing. Thus the final result may be the result of an "
  "averaging with less ",
  "             than 9 points.",
  "",
  "PARAMETER",
  "    nsmooth    INTEGER  Number of times to smooth, default nsmooth=1",
  "    radius     STRING   Search radius, default radius=1deg (units: deg, "
  "rad, km, m)",
  "    maxpoints  INTEGER  Maximum number of points, default "
  "maxpoints=<gridsize>",
  "    form       STRING   Form of the curve, default form=linear",
  "    weight0    FLOAT    Weight at distance 0, default weight0=0.25",
  "    weightR    FLOAT    Weight at the search radius, default weightR=0.25",
};

std::vector<std::string> ReplacevaluesHelp = {
  "NAME",
  "    setvals, setrtoc, setrtoc2 - Replace variable values",
  "",
  "SYNOPSIS",
  "    setvals,oldval,newval[,...]  infile outfile",
  "    setrtoc,rmin,rmax,c  infile outfile",
  "    setrtoc2,rmin,rmax,c,c2  infile outfile",
  "",
  "DESCRIPTION",
  "    This module replaces old variable values with new values, depending on "
  "the operator.",
  "",
  "OPERATORS",
  "    setvals   Set list of old values to new values",
  "              Supply a list of n pairs of old and new values.",
  "    setrtoc   Set range to constant",
  "                       / c      if i(t,x) GE rmin AND i(t,x) LE rmax",
  "              o(t,x) = ",
  "                       \\ i(t,x) if i(t,x) LT rmin AND i(t,x) GT rmax",
  "    setrtoc2  Set range to constant others to constant2",
  "                       / c      if i(t,x) GE rmin AND i(t,x) LE rmax",
  "              o(t,x) = ",
  "                       \\ c2     if i(t,x) LT rmin AND i(t,x) GT rmax",
  "",
  "PARAMETER",
  "    oldval,newval,...  FLOAT   Pairs of old and new values",
  "    rmin               FLOAT   Lower bound",
  "    rmax               FLOAT   Upper bound",
  "    c                  FLOAT   New value - inside range",
  "    c2                 FLOAT   New value - outside range",
};

std::vector<std::string> TimsortHelp = {
  "NAME",
  "    timsort - Timsort",
  "",
  "SYNOPSIS",
  "    timsort  infile outfile",
  "",
  "DESCRIPTION",
  "    Sorts the elements in ascending order over all timesteps for every "
  "field position.",
  "    After sorting it is:",
  "    ",
  "    o(t_1,x) <= o(t_2,x)      forall (t_1<t_2),x",
};

std::vector<std::string> VargenHelp = {
  "NAME",
  "    const, random, topo, for, stdatm - Generate a field",
  "",
  "SYNOPSIS",
  "    const,const,grid  outfile",
  "    random,grid[,seed]  outfile",
  "    topo[,grid]  outfile",
  "    for,start,end[,inc]  outfile",
  "    stdatm,levels  outfile",
  "",
  "DESCRIPTION",
  "    Generates a dataset with one or more fields",
  "",
  "OPERATORS",
  "    const   Create a constant field",
  "            Creates a constant field. All field elements of the grid have "
  "the same value.",
  "    random  Create a field with random numbers",
  "            Creates a field with rectangularly distrubuted random numbers "
  "in the interval [0,1].",
  "    topo    Create a field with topography",
  "            Creates a field with topography data, per default on a global "
  "half degree grid.",
  "    for     Create a time series",
  "            Creates a time series with field size 1 and field elements "
  "beginning with a start value in time step 1",
  "            which is increased from one time step to the next.",
  "    stdatm  Create values for pressure and temperature for hydrostatic "
  "atmosphere",
  "            Creates pressure and temperature values for the given list of "
  "vertical levels.",
  "            The formulars are:",
  "            ",
  "            P(z) = P_0 * exp(-1 * g/R * H/T_0 * log( (exp(z/H)*T_0 + "
  "T_Delta)/(T_0 + T_Delta))",
  "            T(z) = T_0 + T_Delta * exp(-z/H)",
  "            ",
  "            with the following constants",
  "            ",
  "            T_0     = 213 K           Offset to get a surface temperature "
  "of 288K",
  "            T_Delta = 75 K            Temperature lapse rate for 10Km",
  "            P_0     = 1013.25 hPa     Surface pressure",
  "            H       = 10000.0 m       Scale height",
  "            g       = 9.80665 m/s**2  Earth gravity",
  "            R       = 287.05 J/kg*K   Gas constant for air",
  "            ",
  "            This is the solution for the hydrostatic equations and is only "
  "valid for the",
  "            troposphere (constant positive lapse rate). The temperature "
  "increase in the",
  "            stratosphere and other effects of the upper atmosphere are not "
  "taken into",
  "            account.",
  "",
  "PARAMETER",
  "    const   FLOAT   Constant",
  "    seed    INTEGER The seed for a new sequence of pseudo-random numbers "
  "[default: 1]",
  "    grid    STRING  Target grid description file or name",
  "    start   FLOAT   Start value of the loop",
  "    end     FLOAT   End value of the loop",
  "    inc     FLOAT   Increment of the loop [default: 1]",
  "    levels  FLOAT   Target levels in metre above surface",
};

std::vector<std::string> WindTransHelp = {
  "NAME",
  "    uvDestag, rotuvNorth, projuvLatLon - Wind Transformation",
  "",
  "SYNOPSIS",
  "    uvDestag,u,v[,-/+0.5[,-/+0.5]]  infile outfile",
  "    rotuvNorth,u,v  infile outfile",
  "    projuvLatLon,u,v  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains special operators for datsets with wind components "
  "on a rotated lon/lat grid, ",
  "    e.g. data from the regional model HIRLAM or REMO. ",
  "",
  "OPERATORS",
  "    uvDestag      Destaggering of u/v wind components",
  "                  This is a special operator for destaggering of wind "
  "components.",
  "                  If the file contains a grid with temperature (name='t' or "
  "code=11)",
  "                  then grid_temp will be used for destaggered wind.",
  "    rotuvNorth    Rotate u/v wind to North pole.",
  "                  This is an operator for transformation of wind-vectors "
  "from grid-relative to north-pole",
  "                  relative for the whole file. (FAST implementation with "
  "JACOBIANS)",
  "    projuvLatLon  Cylindrical Equidistant projection",
  "                  Thus is an operator for transformation of wind-vectors "
  "from the globe-spherical coordinate system",
  "                  into a flat Cylindrical Equidistant (lat-lon) projection. "
  "(FAST JACOBIAN implementation)",
  "",
  "PARAMETER",
  "    u,v            STRING  Pair of u,v wind components (use variable names "
  "or code numbers)",
  "    -/+0.5,-/+0.5  STRING  Destaggered grid offsets are optional (default "
  "-0.5,-0.5)",
};

std::vector<std::string> RotuvbHelp = {
  "NAME",
  "    rotuvb - Rotation",
  "",
  "SYNOPSIS",
  "    rotuvb,u,v,...  infile outfile",
  "",
  "DESCRIPTION",
  "    This is a special operator for datsets with wind components on a "
  "rotated grid, ",
  "    e.g. data from the regional model REMO. It performs a backward "
  "transformation of ",
  "    velocity components U and V from a rotated spherical system to a "
  "geographical system.",
  "",
  "PARAMETER",
  "    u,v,...  STRING  Pairs of zonal and meridional velocity components (use "
  "variable names or code numbers)",
};

std::vector<std::string> MastrfuHelp = {
  "NAME",
  "    mastrfu - Mass stream function",
  "",
  "SYNOPSIS",
  "    mastrfu  infile outfile",
  "",
  "DESCRIPTION",
  "    This is a special operator for the post processing of the atmospheric "
  "general circulation",
  "    model ECHAM. It computes the mass stream function (code=272). The input "
  "dataset have ",
  "    to be a zonal mean of v-velocity [m/s] (code=132) on pressure levels.",
};

std::vector<std::string> DeriveparHelp = {
  "NAME",
  "    sealevelpressure - Sea level pressure",
  "",
  "SYNOPSIS",
  "    sealevelpressure  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator computes the sea level pressure "
  "(air_pressure_at_sea_level). Required input fields",
  "    are surface_air_pressure, surface_geopotential and air_temperature on "
  "hybrid sigma pressure levels.",
};

std::vector<std::string> AdisitHelp = {
  "NAME",
  "    adisit, adipot - Potential temperature to in-situ temperature and vice "
  "versa",
  "",
  "SYNOPSIS",
  "    adisit[,pressure]  infile outfile",
  "    adipot  infile outfile",
  "",
  "DESCRIPTION",
  "",
  "OPERATORS",
  "    adisit  Potential temperature to in-situ temperature",
  "            This is a special operator for the post processing of the ocean "
  "and sea ice model output.",
  "            It converts potential temperature adiabatically to in-situ "
  "temperature to(t, s, p).",
  "            Required input fields are sea water potential temperature "
  "(name=tho; code=2) and sea water salinity (name=sao; code=5).",
  "            Pressure is calculated from the level information or can be "
  "specified by the optional parameter.",
  "            Output fields are sea water temperature (name=to; code=20) and "
  "sea water salinity (name=s; code=5).",
  "    adipot  In-situ temperature to potential temperature",
  "            This is a special operator for the post processing of the ocean "
  "and sea ice",
  "            model outpu.  It converts in-situ temperature to potential "
  "temperature tho(to,",
  "            s, p).  Required input fields are sea water in-situ temperature "
  "(name=t; code=2) ",
  "            and sea water salinity (name=sao,s; code=5).  Pressure is "
  "calculated",
  "            from the level information or can be specified by the optional "
  "parameter.",
  "            Output fields are sea water temperature (name=tho; code=2) and "
  "sea water",
  "            salinity (name=s; code=5).",
  "",
  "PARAMETER",
  "    pressure  FLOAT   Pressure in bar (constant value assigned to all "
  "levels)",
};

std::vector<std::string> RhopotHelp = {
  "NAME",
  "    rhopot - Calculates potential density",
  "",
  "SYNOPSIS",
  "    rhopot[,pressure]  infile outfile",
  "",
  "DESCRIPTION",
  "    This is a special operator for the post processing of the ocean and sea "
  "ice model MPIOM.",
  "    It calculates the sea water potential density (name=rhopoto; code=18). "
  "Required input fields ",
  "    are sea water in-situ temperature (name=to; code=20) and sea water "
  "salinity (name=sao; code=5).",
  "    Pressure is calculated from the level information or can be specified "
  "by the optional parameter.",
  "",
  "PARAMETER",
  "    pressure  FLOAT   Pressure in bar (constant value assigned to all "
  "levels)",
};

std::vector<std::string> HistogramHelp = {
  "NAME",
  "    histcount, histsum, histmean, histfreq - Histogram",
  "",
  "SYNOPSIS",
  "    <operator>,bounds  infile outfile",
  "",
  "DESCRIPTION",
  "    This module creates bins for a histogram of the input data.",
  "    The bins have to be adjacent and have non-overlapping intervals.",
  "    The user has to define the bounds of the bins. The first value",
  "    is the lower bound and the second value the upper bound of the",
  "    first bin. The bounds of the second bin are defined by the",
  "    second and third value, aso.",
  "    Only 2-dimensional input fields are allowed. The output file ",
  "    contains one vertical level for each of the bins requested.",
  "",
  "OPERATORS",
  "    histcount  Histogram count",
  "               Number of elements in the bin range.",
  "    histsum    Histogram sum",
  "               Sum of elements in the bin range.",
  "    histmean   Histogram mean",
  "               Mean of elements in the bin range.",
  "    histfreq   Histogram frequency",
  "               Relative frequency of elements in the bin range.",
  "",
  "PARAMETER",
  "    bounds  FLOAT  Comma separated list of the bin bounds (-inf and inf "
  "valid)",
};

std::vector<std::string> SethaloHelp = {
  "NAME",
  "    sethalo - Set the left and right bounds of a field",
  "",
  "SYNOPSIS",
  "    sethalo,lhalo,rhalo  infile outfile",
  "",
  "DESCRIPTION",
  "    This operator sets the left and right bounds of the rectangularly "
  "understood fields.",
  "    Positive numbers of the parameter lhalo enlarges the left bound by the "
  "given ",
  "    number of columns from the right bound. The parameter rhalo does the "
  "similar ",
  "    for the right bound. Negative numbers of the parameter lhalo/rhalo can ",
  "    be used to remove the given number of columns of the left and right "
  "bounds.",
  "",
  "PARAMETER",
  "    lhalo  INTEGER  Left halo",
  "    rhalo  INTEGER  Right halo",
};

std::vector<std::string> WctHelp = {
  "NAME",
  "    wct - Windchill temperature",
  "",
  "SYNOPSIS",
  "    wct  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 and infile2 be time series of temperature and wind",
  "    speed records, then a corresponding time series of resulting windchill",
  "    temperatures is written to outfile. The wind chill temperature",
  "    calculation is only valid for a temperature of T <= 33 °C and a wind "
  "speed",
  "    of v >= 1.39 m/s. Whenever these conditions are not satisfied, a "
  "missing",
  "    value is written to outfile. Note that temperature and wind speed "
  "records",
  "    have to be given in units of °C and m/s, respectively.",
};

std::vector<std::string> FdnsHelp = {
  "NAME",
  "    fdns - Frost days where no snow index per time period",
  "",
  "SYNOPSIS",
  "    fdns  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily minimum temperature TN",
  "    and infile2 be a corresponding series of daily surface snow",
  "    amounts. Then the number of days where TN < 0 °C and the surface ",
  "    snow amount is less than 1 cm is counted. The temperature TN",
  "    have to be given in units of Kelvin.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
};

std::vector<std::string> StrwinHelp = {
  "NAME",
  "    strwin - Strong wind days index per time period",
  "",
  "SYNOPSIS",
  "    strwin[,v]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum horizontal wind speed",
  "    VX, then the number of days where VX > v is counted. The horizontal "
  "wind",
  "    speed v is an optional parameter with default v = 10.5 m/s. A further",
  "    output variable is the maximum number of consecutive days with maximum "
  "wind",
  "    speed greater than or equal to v. Note that both VX and v have to be "
  "given in",
  "    units of m/s. Also note that the horizontal wind speed is defined as "
  "the",
  "    square root of the sum of squares of the zonal and meridional wind "
  "speeds.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
  "",
  "PARAMETER",
  "    v  FLOAT   Horizontal wind speed threshold (m/s, default v = 10.5 m/s)",
};

std::vector<std::string> StrbreHelp = {
  "NAME",
  "    strbre - Strong breeze days index per time period",
  "",
  "SYNOPSIS",
  "    strbre  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum horizontal wind speed",
  "    VX, then the number of days where VX is greater than or equal to 10.5 "
  "m/s ",
  "    is counted. A further output variable is the maximum number of "
  "consecutive",
  "    days with maximum wind speed greater than or equal to 10.5 m/s. Note "
  "that",
  "    VX is defined as the square root of the sum of squares of the zonal and",
  "    meridional wind speeds and have to be given in units of m/s.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
};

std::vector<std::string> StrgalHelp = {
  "NAME",
  "    strgal - Strong gale days index per time period",
  "",
  "SYNOPSIS",
  "    strgal  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum horizontal wind speed",
  "    VX, then the number of days where VX is greater than or equal to 20.5 "
  "m/s ",
  "    is counted. A further output variable is the maximum number of "
  "consecutive",
  "    days with maximum wind speed greater than or equal to 20.5 m/s. Note "
  "that",
  "    VX is defined as the square root of the sum of square of the zonal and",
  "    meridional wind speeds and have to be given in units of m/s.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
};

std::vector<std::string> HurrHelp = {
  "NAME",
  "    hurr - Hurricane days index per time period",
  "",
  "SYNOPSIS",
  "    hurr  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum horizontal wind speed",
  "    VX, then the number of days where VX is greater than or equal to 32.5 "
  "m/s",
  "    is counted. A further output variable is the maximum number of "
  "consecutive",
  "    days with maximum wind speed greater than or equal to 32.5 m/s. Note "
  "that",
  "    VX is defined as the square root of the sum of squares of the zonal and",
  "    meridional wind speeds and have to be given in units of m/s.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
};

std::vector<std::string> CMORliteHelp = {
  "NAME",
  "    cmorlite - CMOR lite",
  "",
  "SYNOPSIS",
  "    cmorlite,table[,convert]  infile outfile",
  "",
  "DESCRIPTION",
  "    The CMOR (Climate Model Output Rewriter) library comprises a set of",
  "    functions, that can be used to produce CF-compliant NetCDF files that ",
  "    fulfill the requirements of many of the climate community's standard",
  "    model experiments. These experiments are collectively referred to as",
  "    MIP's. Much of the metadata written to the output files is defined in",
  "    MIP-specific tables, typically made available from each MIP's web site.",
  "    ",
  "    The CDO operator cmorlite process the header and variable section",
  "    of such MIP tables and writes the result with the internal IO library "
  "CDI.",
  "    In addition to the CMOR 2 and 3 table format, the CDO parameter table "
  "format",
  "    is also supported. The following parameter table entries are available:",
  "    ",
  "     Entry           & Type        & Description      ",
  "     name            & WORD        & Name of the variable",
  "     out_name        & WORD        & New name of the variable",
  "     type            & WORD        & Data type (real or double)",
  "     standard_name   & WORD        & As defined in the CF standard name "
  "table",
  "     long_name       & STRING      & Describing the variable",
  "     units           & STRING      & Specifying the units for the variable",
  "     comment         & STRING      & Information concerning the variable",
  "     cell_methods    & STRING      & Information concerning calculation of "
  "means or climatologies",
  "     cell_measures   & STRING      & Indicates the names of the variables "
  "containing cell areas and volumes",
  "     missing_value   & FLOAT       & Specifying how missing data will be "
  "identified",
  "     valid_min       & FLOAT       & Minimum valid value",
  "     valid_max       & FLOAT       & Maximum valid value",
  "     ok_min_mean_abs & FLOAT       & Minimum absolute mean",
  "     ok_max_mean_abs & FLOAT       & Maximum absolute mean",
  "     factor          & FLOAT       & Scale factor",
  "     delete          & INTEGER     & Set to 1 to delete variable",
  "     convert         & INTEGER     & Set to 1 to convert the unit if "
  "necessary",
  "    ",
  "    Most of the above entries are stored as variables attributes, some of "
  "them are handled differently.",
  "    The variable name is used as a search key for the parameter table. "
  "valid_min, valid_max,",
  "    ok_min_mean_abs and ok_max_mean_abs are used to check the range of the "
  "data.",
  "",
  "PARAMETER",
  "    table    STRING   Name of the CMOR table as specified from PCMDI",
  "    convert  STRING   Converts the units if necessary",
};

std::vector<std::string> NCL_windHelp = {
  "NAME",
  "    uv2vr_cfd, uv2dv_cfd - Wind transformation",
  "",
  "SYNOPSIS",
  "    <operator>[,u,v,boundOpt,outMode]  infile outfile",
  "",
  "DESCRIPTION",
  "    This module contains CDO operators with an interface to NCL functions.",
  "    The corresponding NCL functions have the same name. A more detailed "
  "description",
  "    of those NCL function can be found on the NCL homepage "
  "https://www.ncl.ucar.edu.",
  "",
  "OPERATORS",
  "    uv2vr_cfd  U and V wind to relative vorticity",
  "               Computes relative vorticity for a latitude-longitude grid "
  "using centered finite differences.",
  "               The grid need not be global and missing values are allowed.",
  "    uv2dv_cfd  U and V wind to divergence",
  "               Computes divergence for a latitude-longitude grid using "
  "centered finite differences.",
  "               The grid need not be global and missing values are allowed.",
  "",
  "PARAMETER",
  "    u         STRING   Name of variable u (default: u)",
  "    v         STRING   Name of variable v (default: v)",
  "    boundOpt  INTEGER  Boundary condition option (0-3) (default: 0/1 for "
  "cyclic grids)",
  "    outMode   STRING   Output mode new/append (default: new)",
};

std::vector<std::string> CMORHelp = {
  "NAME",
  "    cmor - Climate Model Output Rewriting to produce CMIP-compliant data",
  "",
  "SYNOPSIS",
  "    cmor,MIPtable[,cmor_name=VarList[,key=value[,...]]]  infile",
  "",
  "DESCRIPTION",
  "    ",
  "    ",
  "    The CDO operator cmor converts an infile into a CMIP-compliant format",
  "    by using the CMOR library. Each output file contains a single output "
  "variable.",
  "    The name of the output files are generated by CMOR based on the DRS "
  "(Data reference",
  "    Syntax) of the project. CMOR checks and applies the information "
  "delivered",
  "    through the project dependend MIPtable on the infile. Additional "
  "information",
  "    which is required for the conversion can be configured via keyvalues as "
  "optional parameters.",
  "    ",
  "    By specifying a variable selector keyvalue, e.g. cmor_name=tas, the "
  "user can",
  "    pre-select a subset of infile variables. If name or code is specified, "
  "a",
  "    corresponding cmor_name which can also be found in the MIPtable is also",
  "    required to map the infile variable to the CMOR-variable. For mapping "
  "more",
  "    variables at the operator call, one can specify a mapping table via "
  "keyword mapping_table.",
  "    ",
  "    Global attributes must be collected in info files and can be specified "
  "via keyword",
  "    info. All required and optional global attributes as well as "
  "information",
  "    about table file formats are given in the 'cdo cmor manual'.",
  "    ",
  "    If questions remain, do not hesitate to ask and send an email to "
  "wachsmannATdkrz.de.",
  "    ",
  "",
  "PARAMETER",
  "    MIPtable                   STRING    Name of the MIP table as used by "
  "CMOR.",
  "                               "
  "----------------------------------------------------------------------------"
  "----------------",
  "    cmor_name           | cn   STRING    Variable selector and specified in "
  "the MIP table.",
  "                                         Comma separated list of CMOR "
  "variable names.",
  "                                         Default is to process all "
  "variables.",
  "    name                | n    STRING    Variable selector.",
  "                                         Name of a selected @file{infile} "
  "variable.",
  "    code                | c    INTEGER   Variable selector. ",
  "                                         Three digits (GRIB) Code of a "
  "selected @file{infile} variable.",
  "                               "
  "----------------------------------------------------------------------------"
  "----------------",
  "    info                | i    STRING    Preprozessing.",
  "                                         List of filenames containing "
  "global attributes and information.",
  "                                         Default: CWD/.cdocmorinfo",
  "    grid_info           | gi   STRING    Preprozessing.",
  "                                         NetCDF file with model grid "
  "description.",
  "    mapping_table       | mt   STRING    Preprozessing.",
  "                                         Fortran Namelist containing "
  "variable information for e.g. renaming.",
  "                               "
  "----------------------------------------------------------------------------"
  "----------------",
  "    drs                 | d    CHARACTER Output control.",
  "                                         Do(=y, default) or do not(=n) move "
  "output into the project DRS structure.",
  "    drs_root            | dr   STRING    Output control. CMOR output root "
  "directory.",
  "                                         Default: CWD.",
  "    output_mode         | om   CHARACTER Output control.",
  "                                         Either 'r' for replace (default) "
  "or 'a' for append mode.",
  "    last_chunk          | lc   STRING    Output control. Filename of chunk "
  "to which shall be appended.  ",
  "    maximum_size        | ms   INTEGER   Output control. Limit of output "
  "file sie in GigaByte.",
  "                               "
  "----------------------------------------------------------------------------"
  "----------------",
  "    required_time_units | rtu  STRING    Temporal description.",
  "                                         Time axis reference date specified "
  "by the experiment.",
  "                                         Format: 'days since YYYY-day-month "
  "hh:mm:ss'.",
  "    cell_methods        | cm   CHARACTER Temporal description.",
  "                                         Cell_methods of time axis.",
  "                                         Value is one of 'm' (default)  , "
  "'p', 'c', 'n'.",
  "                               "
  "----------------------------------------------------------------------------"
  "----------------",
  "    units               | u    STRING    Variable attrbiute. Units of the "
  "variable.",
  "                                         Must be known by library UDunits.",
  "    variable_comment    | vc   STRING    Variable attrbiute. Variable "
  "comment.",
  "    positive            | p    CHARACTER Variable attrbiute.",
  "                                         Positive flux direction, either "
  "'u' for upward or 'd' for downward.",
  "                               "
  "----------------------------------------------------------------------------"
  "----------------",
  "    scalar_z_coordinate | szc  STRING    Scalar z-coordinate (=name_value).",
  "    character_axis      | ca   STRING    CMOR name of a character axis.",
  "                                         Valid axes are: basin, vegtype or "
  "oline. ",
};

std::vector<std::string> MagplotHelp = {
  "NAME",
  "    contour, shaded, grfill - Lat/Lon plot",
  "",
  "SYNOPSIS",
  "    <operator>,params  infile obase",
  "",
  "DESCRIPTION",
  "    The operators in this module generates 2D Lon/Lat plots.",
  "    The data for the plot is read from infile.",
  "    Only data on rectilinear Lon/Lat grids are supported.",
  "    The output file will be named <obase>_<param>.<device> where param is "
  "the parameter name and",
  "    device is the device name. The default output file format is "
  "postscript,",
  "    this can be changed with the device parameter.",
  "    The type of the plot depends on the choosen operator.",
  "    ",
  "    Here is a list of all common plot parameters:",
  "    ",
  "     Keyname     & Type    & Description      ",
  "     device      & STRING  & Output device (ps, eps, pdf, png, gif, "
  "gif_animation, jpeg, svg, kml)",
  "     projection  & STRING  & Projection (cylindrical, polar_stereographic, "
  "robinson, mercator)",
  "     style       & STRING  & Contour line style (solid, dash, dot, "
  "chain_dash, chain_dot)",
  "     min         & FLOAT   & Minimum value",
  "     max         & FLOAT   & Maximum value",
  "     lon_max     & FLOAT   & Maximum longitude of the image",
  "     lon_min     & FLOAT   & Minimum longitude of the image",
  "     lat_max     & FLOAT   & Maximum latitude of the image",
  "     lat_min     & FLOAT   & Minimum latitude of the image",
  "     count       & INTEGER & Number of Contour levels / Colour bands  ",
  "     interval    & FLOAT   & Interval in data units between two bands lines",
  "     list        & INTEGER & List of levels to be plotted",
  "     RGB         & STRING  & TRUE or FALSE, to  indicate, if the input "
  "colour is in RGB format",
  "     step_freq   & INTEGER & Frequency of time steps to be considered for "
  "making the animation",
  "                 &         & (device=gif_animation). Default value is \"1\" "
  "(all time steps).",
  "                 &         & Will be ignored if input file has multiple "
  "variables.",
  "     file_split  & STRING  & TRUE or FALSE, to split the output file for "
  "each variable, if input has",
  "                 &         & multiple variables. Default value is "
  "\"FALSE\". Valid only for \"PS\" format.",
  "",
  "OPERATORS",
  "    contour  Contour plot",
  "             The operator contour generates the discrete contour lines of "
  "the input field values.",
  "             The following additional parameters are valid for contour "
  "operator,",
  "             module in addition to the common plot parameters:",
  "             ",
  "              Keyname      & Type    & Description      ",
  "              colour       & STRING  & Colour for drawing the contours",
  "              thickness    & FLOAT   & Thickness of the contour line",
  "              style        & STRING  & Line Style can be \"SOLID\", "
  "\"DASH\", \"DOT\", \"CHAIN_DASH\",",
  "                           &         & \"CHAIN_DOT\"",
  "    shaded   Shaded contour plot",
  "             The operator shaded generates the filled contours of the given "
  "input field values.",
  "             The following additional parameters are valid for shaded "
  "contour and gridfill operator,",
  "             in addition to the common plot parameters.",
  "             ",
  "              Keyname      & Type    & Description      ",
  "              colour_min   & STRING  & Colour for the Minimum colour band",
  "              colour_max   & STRING  & Colour for the Minimum colour band",
  "              colour_triad & STRING  & Direction of colour sequencing for "
  "shading \"CW\" or \"ACW\",",
  "                           &         & to denote \"clockwise\" and "
  "\"anticlockwise\" respectively.",
  "                           &         & To be used in conjunction with "
  "\"colour_min\", \"colour_max\"",
  "                           &         & options. Default is \"ACW\"",
  "              colour_table & STRING  & File with user specified colours "
  "with the format as",
  "             ",
  "             Example file for 6 colours in RGB format:",
  "             	6",
  "             	RGB(0.0;0.0;1.0)",
  "             	RGB(0.0;0.0;0.5)",
  "             	RGB(0.0;0.5;0.5)",
  "             	RGB(0.0;1.0;0.0)",
  "             	RGB(0.5;0.5;0.0)",
  "             	RGB(1.0;0.0;0.0)",
  "             ",
  "    grfill   Shaded gridfill plot",
  "             The operator grfill is similar to satellite imaging and shades "
  "each cell (pixel) according",
  "             to the value of the field at that cell.",
  "",
  "PARAMETER",
  "    params  STRING   Comma separated list of plot parameters",
  "",
  "NOTE",
  "    All colour parameter can be either standard name or in RGB format.",
  "    The valid standard name strings for \"colour\" are:",
  "    ",
  "    \"red\", \"green\", \"blue\", \"yellow\", \"cyan\", \"magenta\", "
  "\"black\", \"avocado\", \"beige\",",
  "    \"brick\", \"brown\", \"burgundy\", \"charcoal\", \"chestnut\", "
  "\"coral\", \"cream\", \"evergreen\",",
  "    \"gold\", \"grey\", \"khaki\", \"kellygreen\", \"lavender\", "
  "\"mustard\", \"navy\", \"ochre\",",
  "    \"olive\", \"peach\", \"pink\", \"rose\", \"rust\", \"sky\", \"tan\", "
  "\"tangerine\", \"turquoise\",",
  "    \"violet\", \"reddishpurple\", \"purplered\", \"purplishred\", "
  "\"orangishred\", \"redorange\",",
  "    \"reddishorange\", \"orange\", \"yellowishorange\", \"orangeyellow\", "
  "\"orangishyellow\",",
  "    \"greenishyellow\", \"yellowgreen\", \"yellowishgreen\", "
  "\"bluishgreen\", \"bluegreen\",",
  "    \"greenishblue\", \"purplishblue\", \"bluepurple\", \"bluishpurple\", "
  "\"purple\", \"white\"",
};

std::vector<std::string> MagvectorHelp = {
  "NAME",
  "    vector - Lat/Lon vector plot",
  "",
  "SYNOPSIS",
  "    vector,params  infile obase",
  "",
  "DESCRIPTION",
  "    This operator generates 2D Lon/Lat vector plots.",
  "    The data for the plot is read from infile. The input is expected to "
  "contain two velocity",
  "    components. Only data on rectilinear Lon/Lat grids are supported.",
  "    The output file will be named <obase>.<device> where device is the "
  "device name. ",
  "    The default output file format is postscript, this can be changed with "
  "the device parameter.",
  "    ",
  "    Here is a list of all vector plot parameters:",
  "    ",
  "     Keyname     & Type    & Description      ",
  "     device      & STRING  & Output device (ps, eps, pdf, png, gif, "
  "gif_animation, jpeg, svg, kml)",
  "     projection  & STRING  & Projection (cylindrical, polar_stereographic, "
  "robinson, mercator)",
  "     thin_fac    & FLOAT   & Controls the actual number of wind arrows or "
  "flags plotted (default 2).",
  "     unit_vec    & FLOAT   & Wind speed in m/s represented by a unit vector "
  "(1.0cm)",
  "     step_freq   & INTEGER & Frequency of time steps to be considered for "
  "making the animation",
  "                 &         & (device=gif_animation). Default value is \"1\" "
  "(all time steps).",
  "                 &         & Will be ignored if input file has multiple "
  "variables.",
  "",
  "PARAMETER",
  "    params  STRING   Comma separated list of plot parameters",
};

std::vector<std::string> MaggraphHelp = {
  "NAME",
  "    graph - Line graph plot",
  "",
  "SYNOPSIS",
  "    graph,params  infiles outfile",
  "",
  "DESCRIPTION",
  "    This operator generates line graph plots.",
  "    The data for the plot is read from infiles. The result is written to "
  "outfile.",
  "    The default output file format is postscript, this can be changed with "
  "the device parameter.",
  "    ",
  "    Here is a list of all graph plot parameters:",
  "    ",
  "     Keyname    & Type    & Description      ",
  "     device     & STRING  & Output device (ps, eps, pdf, png, gif, "
  "gif_animation, jpeg, svg, kml)",
  "     ymin       & FLOAT   & Minimum value of the y-axis data ",
  "     ymax       & FLOAT   & Maximum value of the y-axis data ",
  "     linewidth  & INT     & Linewidth (default 8)",
  "     stat       & STRING  & \"TRUE\" or \"FALSE\", to switch on the mean "
  "computation. Default is \"FALSE\".",
  "                &         & Will be overridden to \"FALSE\", if input files "
  "have unequal number of time",
  "                &         & steps or different start/end times. ",
  "     sigma      & FLOAT   & Standard deviation value for generating shaded "
  "back ground around the mean value.",
  "                &         & To be used in conjunction with 'stat=\"TRUE\"' ",
  "     obsv       & STRING  & To indicate if the input files have an "
  "observation data, by setting to \"TRUE\".",
  "                &         & Default value is \"FALSE\". The observation "
  "data should be the first file in the",
  "                &         & input file list. The observation data is always "
  "plotted in black colour. ",
  "",
  "PARAMETER",
  "    params  STRING   Comma separated list of plot parameters",
};

std::vector<std::string> EcaCddHelp = {
  "NAME",
  "    eca_cdd - Consecutive dry days index per time period",
  "",
  "SYNOPSIS",
  "    eca_cdd[,R[,N]]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily precipitation amount RR, then "
  "the largest number ",
  "    of consecutive days where RR is less than R is counted. R is an "
  "optional parameter with ",
  "    default R = 1 mm. A further output variable is the number of dry "
  "periods of more than N days.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    R  FLOAT    Precipitation threshold (unit: mm; default: R = 1 mm)",
  "    N  INTEGER  Minimum number of days exceeded (default: N = 5)",
};

std::vector<std::string> EcaCfdHelp = {
  "NAME",
  "    eca_cfd - Consecutive frost days index per time period",
  "",
  "SYNOPSIS",
  "    eca_cfd[,N]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily minimum temperature TN, then "
  "the largest number of",
  "    consecutive days where TN < 0 °C is counted. Note that TN have to be "
  "given in units of Kelvin.",
  "    A further output variable is the number of frost periods of more than N "
  "days.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    N  INTEGER  Minimum number of days exceeded (default: N = 5)",
};

std::vector<std::string> EcaCsuHelp = {
  "NAME",
  "    eca_csu - Consecutive summer days index per time period",
  "",
  "SYNOPSIS",
  "    eca_csu[,T[,N]]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum temperature TX, then "
  "the largest number of consecutive",
  "    days where TX > T is counted. The number T is an optional parameter "
  "with default T = 25°C.",
  "    Note that TN have to be given in units of Kelvin, whereas T have to be "
  "given in degrees Celsius.",
  "    A further output variable is the number of summer periods of more than "
  "N days.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    T  FLOAT    Temperature threshold (unit: °C; default: T = 25°C)",
  "    N  INTEGER  Minimum number of days exceeded (default: N = 5)",
};

std::vector<std::string> EcaCwdHelp = {
  "NAME",
  "    eca_cwd - Consecutive wet days index per time period",
  "",
  "SYNOPSIS",
  "    eca_cwd[,R[,N]]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily precipitation amount RR, then "
  "the largest number ",
  "    of consecutive days where RR is at least R is counted. R is an optional "
  "parameter with ",
  "    default R = 1 mm. A further output variable is the number of wet "
  "periods of more than N days.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    R  FLOAT    Precipitation threshold (unit: mm; default: R = 1 mm)",
  "    N  INTEGER  Minimum number of days exceeded (default: N = 5)",
};

std::vector<std::string> EcaCwdiHelp = {
  "NAME",
  "    eca_cwdi - Cold wave duration index wrt mean of reference period",
  "",
  "SYNOPSIS",
  "    eca_cwdi[,nday[,T]]  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily minimum temperature TN, and "
  "let infile2 be the mean ",
  "    TNnorm of daily minimum temperatures for any period used as reference. "
  "Then counted is the number of days",
  "    where, in intervals of at least nday consecutive days, TN < TNnorm - T.",
  "    The numbers nday and T are optional parameters with default nday = 6 "
  "and T = 5°C. ",
  "    A further output variable is the number of cold waves longer than or "
  "equal to nday days.",
  "    TNnorm is calculated as the mean of minimum temperatures of a five day "
  "window centred on each calendar day ",
  "    of a given climate reference period. Note that both TN and TNnorm have "
  "to be given in the same units.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
  "",
  "PARAMETER",
  "    nday  INTEGER  Number of consecutive days (default: nday = 6)",
  "    T     FLOAT    Temperature offset (unit: °C; default: T = 5°C)",
};

std::vector<std::string> EcaCwfiHelp = {
  "NAME",
  "    eca_cwfi - Cold-spell days index wrt 10th percentile of reference "
  "period",
  "",
  "SYNOPSIS",
  "    eca_cwfi[,nday]  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily mean temperature TG, and "
  "infile2 be the 10th",
  "    percentile TGn10 of daily mean temperatures for any period used as "
  "reference. ",
  "    Then counted is the number of days where, in intervals of at least nday "
  "consecutive days,",
  "    TG < TGn10. The number nday is an optional parameter with default nday "
  "= 6.",
  "    A further output variable is the number of cold-spell periods longer "
  "than or equal to nday days.",
  "    TGn10 is calculated as the 10th percentile of daily mean temperatures "
  "of a five day window ",
  "    centred on each calendar day of a given climate reference period. Note "
  "that both TG and TGn10 ",
  "    have to be given in the same units. The date information of a timestep "
  "in outfile is the date of",
  "    the last contributing timestep in infile1.",
  "",
  "PARAMETER",
  "    nday  INTEGER  Number of consecutive days (default: nday = 6)",
};

std::vector<std::string> EcaEtrHelp = {
  "NAME",
  "    eca_etr - Intra-period extreme temperature range",
  "",
  "SYNOPSIS",
  "    eca_etr  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 and infile2 be time series of thr maximum and minimum",
  "    temperature TX and TN, respectively. Then the extreme temperature",
  "    range is the difference of the maximum of TX and the minimum of TN.",
  "    Note that TX and TN have to be given in the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timesteps in infile1 and infile2.",
};

std::vector<std::string> EcaFdHelp = {
  "NAME",
  "    eca_fd - Frost days index per time period",
  "",
  "SYNOPSIS",
  "    eca_fd  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily minimum temperature TN,",
  "    then the number of days where TN < 0 °C is counted. Note",
  "    that TN have to be given in units of Kelvin.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
};

std::vector<std::string> EcaGslHelp = {
  "NAME",
  "    eca_gsl - Thermal Growing season length index",
  "",
  "SYNOPSIS",
  "    eca_gsl[,nday[,T[,fland]]]  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily mean temperature TG, and "
  "infile2 be a land-water mask.",
  "    Within a period of 12 months, the thermal growing season length is "
  "officially defined as the number of days between:",
  "    - first occurrence of at least nday consecutive days with TG > T",
  "    - first occurrence of at least nday consecutive days with TG < T within "
  "the last 6 months",
  "    On northern hemisphere, this period corresponds with the regular year, "
  "whereas on southern hemisphere, it starts ",
  "    at July 1st. Please note, that this definition may lead to weird "
  "results concerning values TG = T: ",
  "    In the first half of the period, these days do not contribute to the "
  "gsl, but they do within the second half.",
  "    Moreover this definition could lead to discontinuous values in "
  "equatorial regions.",
  "    ",
  "    The numbers nday and T are optional parameter with default nday = 6 and "
  "T = 5°C. ",
  "    The number fland is an optional parameter with default value fland = "
  "0.5 and denotes the fraction of ",
  "    a grid point that have to be covered by land in order to be included in "
  "the calculation. A further output variable ",
  "    is the start day of year of the growing season. Note that TG have to be "
  "given in units of Kelvin, whereas T ",
  "    have to be given in degrees Celsius.",
  "    ",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    nday   INTEGER  Number of consecutive days (default: nday = 6)",
  "    T      FLOAT    Temperature threshold (unit: °C; default: T = 5°C)",
  "    fland  FLOAT    Land fraction threshold (default: fland = 0.5)",
};

std::vector<std::string> EcaHdHelp = {
  "NAME",
  "    eca_hd - Heating degree days per time period",
  "",
  "SYNOPSIS",
  "    eca_hd[,T1[,T2]]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily mean temperature TG, then the "
  "heating degree days ",
  "    are defined as the sum of T1 - TG, where only values TG < T2 are "
  "considered. ",
  "    If T1 and T2 are omitted, a temperature of 17°C is used for both "
  "parameters. ",
  "    If only T1 is given, T2 is set to T1. Note that TG have to be given in "
  "units ",
  "    of kelvin, whereas T1 and T2 have to be given in degrees Celsius.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    T1  FLOAT   Temperature limit (unit: °C; default: T1 = 17°C)",
  "    T2  FLOAT   Temperature limit (unit: °C; default: T2 = T1)",
};

std::vector<std::string> EcaHwdiHelp = {
  "NAME",
  "    eca_hwdi - Heat wave duration index wrt mean of reference period",
  "",
  "SYNOPSIS",
  "    eca_hwdi[,nday[,T]]  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily maximum temperature TX, and "
  "let infile2 be the mean ",
  "    TXnorm of daily maximum temperatures for any period used as reference. "
  "Then counted is the number of days",
  "    where, in intervals of at least nday consecutive days, TX > TXnorm + T.",
  "    The numbers nday and T are optional parameters with default nday = 6 "
  "and T = 5°C. ",
  "    A further output variable is the number of heat waves longer than or "
  "equal to nday days. ",
  "    TXnorm is calculated as the mean of maximum temperatures of a five day "
  "window centred on each calendar day",
  "    of a given climate reference period. Note that both TX and TXnorm have "
  "to be given in the same units.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
  "",
  "PARAMETER",
  "    nday  INTEGER  Number of consecutive days (default: nday = 6)",
  "    T     FLOAT    Temperature offset (unit: °C; default: T = 5°C)",
};

std::vector<std::string> EcaHwfiHelp = {
  "NAME",
  "    eca_hwfi - Warm spell days index wrt 90th percentile of reference "
  "period",
  "",
  "SYNOPSIS",
  "    eca_hwfi[,nday]  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily mean temperature TG, and ",
  "    infile2 be the 90th percentile TGn90 of daily mean temperatures",
  "    for any period used as reference. Then counted is the number of days",
  "    where, in intervals of at least nday consecutive days, TG > TGn90. The",
  "    number nday is an optional parameter with default nday = 6. A further",
  "    output variable is the number of warm-spell periods longer than or",
  "    equal to nday days. ",
  "    TGn90 is calculated as the 90th percentile of daily mean temperatures "
  "of a five ",
  "    day window centred on each calendar day of a given climate reference "
  "period.",
  "    Note that both TG and TGn90 have to be given in the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile1.",
  "",
  "PARAMETER",
  "    nday  INTEGER  Number of consecutive days (default: nday = 6)",
};

std::vector<std::string> EcaIdHelp = {
  "NAME",
  "    eca_id - Ice days index per time period",
  "",
  "SYNOPSIS",
  "    eca_id  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum temperature TX,",
  "    then the number of days where TX < 0 °C is counted. Note",
  "    that TX have to be given in units of Kelvin.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
};

std::vector<std::string> EcaR75pHelp = {
  "NAME",
  "    eca_r75p - Moderate wet days wrt 75th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_r75p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 75th percentile RRn75 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the percentage of wet days with RR > RRn75 is "
  "calculated. ",
  "    RRn75 is calculated as the 75th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,75.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR75ptotHelp = {
  "NAME",
  "    eca_r75ptot - Precipitation percent due to R75p days",
  "",
  "SYNOPSIS",
  "    eca_r75ptot  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 75th percentile RRn75 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the ratio of the precipitation sum at wet days "
  "with RR > RRn75 to the total ",
  "    precipitation sum is calculated. ",
  "    RRn75 is calculated as the 75th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,75.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR90pHelp = {
  "NAME",
  "    eca_r90p - Wet days wrt 90th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_r90p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 90th percentile RRn90 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the percentage of wet days with RR > RRn90 is "
  "calculated. ",
  "    RRn90 is calculated as the 90th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,90.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR90ptotHelp = {
  "NAME",
  "    eca_r90ptot - Precipitation percent due to R90p days",
  "",
  "SYNOPSIS",
  "    eca_r90ptot  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 90th percentile RRn90 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the ratio of the precipitation sum at wet days "
  "with RR > RRn90 to the total ",
  "    precipitation sum is calculated. ",
  "    RRn90 is calculated as the 90th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,90.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR95pHelp = {
  "NAME",
  "    eca_r95p - Very wet days wrt 95th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_r95p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 95th percentile RRn95 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the percentage of wet days with RR > RRn95 is "
  "calculated. ",
  "    RRn95 is calculated as the 95th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,95.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR95ptotHelp = {
  "NAME",
  "    eca_r95ptot - Precipitation percent due to R95p days",
  "",
  "SYNOPSIS",
  "    eca_r95ptot  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 95th percentile RRn95 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the ratio of the precipitation sum at wet days "
  "with RR > RRn95 to the total ",
  "    precipitation sum is calculated. ",
  "    RRn95 is calculated as the 95th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,95.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR99pHelp = {
  "NAME",
  "    eca_r99p - Extremely wet days wrt 99th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_r99p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 99th percentile RRn99 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the percentage of wet days with RR > RRn99 is "
  "calculated. ",
  "    RRn99 is calculated as the 99th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,99.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaR99ptotHelp = {
  "NAME",
  "    eca_r99ptot - Precipitation percent due to R99p days",
  "",
  "SYNOPSIS",
  "    eca_r99ptot  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series RR of the daily precipitation amount at "
  "wet days (precipitation >= 1 mm)",
  "    and infile2 be the 99th percentile RRn99 of the daily precipitation "
  "amount at wet days for any period ",
  "    used as reference. Then the ratio of the precipitation sum at wet days "
  "with RR > RRn99 to the total ",
  "    precipitation sum is calculated. ",
  "    RRn99 is calculated as the 99th percentile of all wet days of a given "
  "climate reference period.",
  "    Usually infile2 is generated by the operator ydaypctl,99.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaPdHelp = {
  "NAME",
  "    eca_pd, eca_r10mm, eca_r20mm - Precipitation days index per time period",
  "",
  "SYNOPSIS",
  "    eca_pd,x  infile outfile",
  "    eca_r10mm  infile outfile",
  "    eca_r20mm  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily precipitation amount RR in "
  "[mm] (or alternatively in [kg m-2]),",
  "    then the number of days where RR is at least x mm is counted. ",
  "    eca_r10mm and eca_r20mm are specific ECA operators with a daily "
  "precipitation amount of 10 and 20 mm respectively.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "OPERATORS",
  "    eca_pd     Precipitation days index per time period",
  "               Generic ECA operator with daily precipitation sum exceeding "
  "x mm.",
  "    eca_r10mm  Heavy precipitation days index per time period",
  "               Specific ECA operator with daily precipitation sum exceeding "
  "10 mm.",
  "    eca_r20mm  Very heavy precipitation days index per time period",
  "               Specific ECA operator with daily precipitation sum exceeding "
  "20 mm.",
  "",
  "PARAMETER",
  "    x  FLOAT   Daily precipitation amount threshold in [mm]",
  "",
  "NOTE",
  "    Precipitation rates in [mm/s] have to be converted to precipitation "
  "amounts (multiply with 86400 s).",
  "    Apart from metadata information the result of eca_pd,1 and eca_rr1 is "
  "the same.",
};

std::vector<std::string> EcaRr1Help = {
  "NAME",
  "    eca_rr1 - Wet days index per time period",
  "",
  "SYNOPSIS",
  "    eca_rr1[,R]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily precipitation amount RR in "
  "[mm] (or alternatively in [kg m-2]), then",
  "    the number of days where RR is at least R is counted. R is an optional "
  "parameter with default R = 1 mm. ",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    R  FLOAT   Precipitation threshold (unit: mm; default: R = 1 mm)",
};

std::vector<std::string> EcaRx1dayHelp = {
  "NAME",
  "    eca_rx1day - Highest one day precipitation amount per time period",
  "",
  "SYNOPSIS",
  "    eca_rx1day[,mode]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily precipitation amount RR,",
  "    then the maximum of RR is written to outfile. If the optional",
  "    parameter mode is set to 'm' the maximum daily precipitation",
  "    amounts are determined for each month. ",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile.",
  "",
  "PARAMETER",
  "    mode  STRING   Operation mode (optional). If mode = 'm' then maximum "
  "daily precipitation amounts are determined for each month",
};

std::vector<std::string> EcaRx5dayHelp = {
  "NAME",
  "    eca_rx5day - Highest five-day precipitation amount per time period",
  "",
  "SYNOPSIS",
  "    eca_rx5day[,x]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of 5-day precipitation totals RR, then the "
  "maximum of RR is written to outfile. ",
  "    A further output variable is the number of 5 day period with "
  "precipitation totals greater than x mm, where x ",
  "    is an optional parameter with default x = 50 mm.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    x  FLOAT   Precipitation threshold (unit: mm; default: x = 50 mm)",
};

std::vector<std::string> EcaSdiiHelp = {
  "NAME",
  "    eca_sdii - Simple daily intensity index per time period",
  "",
  "SYNOPSIS",
  "    eca_sdii[,R]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily precipitation amount RR, then "
  "the mean precipitation amount at ",
  "    wet days (RR > R) is written to outfile. R is an optional parameter "
  "with default R = 1 mm.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    R  FLOAT   Precipitation threshold (unit: mm; default: R = 1 mm)",
};

std::vector<std::string> EcaSuHelp = {
  "NAME",
  "    eca_su - Summer days index per time period",
  "",
  "SYNOPSIS",
  "    eca_su[,T]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily maximum temperature TX, then "
  "the number of days where ",
  "    TX > T is counted. The number T is an optional parameter with default T "
  "= 25°C. ",
  "    Note that TX have to be given in units of Kelvin, whereas T have to be "
  "given in degrees Celsius.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    T  FLOAT   Temperature threshold (unit: °C; default: T = 25°C)",
};

std::vector<std::string> EcaTg10pHelp = {
  "NAME",
  "    eca_tg10p - Cold days percent wrt 10th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_tg10p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily mean temperature TG, and",
  "    infile2 be the 10th percentile TGn10 of daily mean temperatures",
  "    for any period used as reference. Then the percentage of time where ",
  "    TG < TGn10 is calculated.",
  "    TGn10 is calculated as the 10th percentile of daily mean temperatures "
  "of a five ",
  "    day window centred on each calendar day of a given climate reference "
  "period.",
  "    Note that both TG and TGn10 have to be given in the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile1.",
};

std::vector<std::string> EcaTg90pHelp = {
  "NAME",
  "    eca_tg90p - Warm days percent wrt 90th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_tg90p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily mean temperature TG, and",
  "    infile2 be the 90th percentile TGn90 of daily mean temperatures",
  "    for any period used as reference. Then the percentage of time where TG "
  "> TGn90 ",
  "    is calculated. ",
  "    TGn90 is calculated as the 90th percentile of daily mean temperatures "
  "of a five ",
  "    day window centred on each calendar day of a given climate reference "
  "period.",
  "    Note that both TG and TGn90 have to be given in the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile1.",
};

std::vector<std::string> EcaTn10pHelp = {
  "NAME",
  "    eca_tn10p - Cold nights percent wrt 10th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_tn10p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time serie of the daily minimum temperature TN, and",
  "    infile2 be the 10th percentile TNn10 of daily minimum temperatures",
  "    for any period used as reference. Then the percentage of time where TN "
  "< TNn10 ",
  "    is calculated.",
  "    TNn10 is calculated as the 10th percentile of daily minimum "
  "temperatures of a five ",
  "    day window centred on each calendar day of a given climate reference "
  "period.",
  "    Note that both TN and TNn10 have to be given in the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile1.",
};

std::vector<std::string> EcaTn90pHelp = {
  "NAME",
  "    eca_tn90p - Warm nights percent wrt 90th percentile of reference period",
  "",
  "SYNOPSIS",
  "    eca_tn90p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily minimum temperature TN, and "
  "infile2 be the ",
  "    90th percentile TNn90 of daily minimum temperatures for any period used "
  "as reference. ",
  "    Then the percentage of time where TN > TNn90 is calculated. TNn90 is "
  "calculated as the 90th percentile",
  "    of daily minimum temperatures of a five day window centred on each "
  "calendar day of a given climate",
  "    reference period. Note that both TN and TNn90 have to be given in the "
  "same units.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile1.",
};

std::vector<std::string> EcaTrHelp = {
  "NAME",
  "    eca_tr - Tropical nights index per time period",
  "",
  "SYNOPSIS",
  "    eca_tr[,T]  infile outfile",
  "",
  "DESCRIPTION",
  "    Let infile be a time series of the daily minimum temperature TN, then "
  "the number of days where ",
  "    TN > T is counted. The number T is an optional parameter with default T "
  "= 20°C. ",
  "    Note that TN have to be given in units of Kelvin, whereas T have to be "
  "given in degrees Celsius.",
  "    The date information of a timestep in outfile is the date of the last "
  "contributing timestep in infile.",
  "",
  "PARAMETER",
  "    T  FLOAT   Temperature threshold (unit: °C; default: T = 20°C)",
};

std::vector<std::string> EcaTx10pHelp = {
  "NAME",
  "    eca_tx10p - Very cold days percent wrt 10th percentile of reference "
  "period",
  "",
  "SYNOPSIS",
  "    eca_tx10p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily maximum temperature TX, and",
  "    infile2 be the 10th percentile TXn10 of daily maximum temperatures",
  "    for any period used as reference. Then the percentage of time where TX "
  "< TXn10.",
  "    is calculated.",
  "    TXn10 is calculated as the 10th percentile of daily maximum "
  "temperatures of a five ",
  "    day window centred on each calendar day of a given climate reference "
  "period.",
  "    Note that both TX and TXn10 have to be givenin the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile1.",
};

std::vector<std::string> EcaTx90pHelp = {
  "NAME",
  "    eca_tx90p - Very warm days percent wrt 90th percentile of reference "
  "period",
  "",
  "SYNOPSIS",
  "    eca_tx90p  infile1 infile2 outfile",
  "",
  "DESCRIPTION",
  "    Let infile1 be a time series of the daily maximum temperature TX, and",
  "    infile2 be the 90th percentile TXn90 of daily maximum temperatures",
  "    for any period used as reference. Then the percentage of time where TX "
  "> TXn90.",
  "    is calculated.",
  "    TXn90 is calculated as the 90th percentile of daily maximum "
  "temperatures of a five ",
  "    day window centred on each calendar day of a given climate reference "
  "period.",
  "    Note that both TX and TXn90 have to be given in the same units.",
  "    The date information of a timestep in outfile is the date of",
  "    the last contributing timestep in infile1.",
};
