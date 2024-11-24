#! julia

push!(LOAD_PATH, pwd()*"/src/")

using Printf
using CSV, DataFrames
using PrettyTables


function  do_it( ARGS )
    df = CSV.read( ARGS[1], DataFrame, types=String );

    println( df );

    ha = LatexHighlighter((d,i,j)->i % 2 != 0,
                          ["cellcolor{lightgray}","texttt"])
    hb = LatexHighlighter((d,i,j)->i % 2 == 0,
                          ["texttt"])

    iox = open(ARGS[2], "w");
    pretty_table( iox,df, header = names( df ),
                  backend = Val(:latex), highlighters = (ha, hb));

    #show( iox, "text/latex", df );
    close( iox );
end

num_args = length( ARGS );
if   num_args != 2
    println( "\n\tcsv2ltex [input.csv] [output.tex]\n\n\n" );
    exit( -1 );
end

do_it( ARGS );
