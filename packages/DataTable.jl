""" # Julia Data Table Package
    Author: Ken Sible | Last Modified: March 8, 2019
    Exported Function(s)  : append!, insert!, remove!
    Exported Structure(s) : Table
"""
module DataTable
import Base: show, append!, insert!
export Table, append!, insert!, remove!

mutable struct Table
    name::String
    clms::Vector{String}
    rows::Matrix
    clm_widths::Vector{Integer}
    max_width::Integer

    function Table(name, clms, rows=Matrix(undef, 0, 0); max_width=20)
        clm_widths = zeros(length(clms))
        for i = 1:length(clms)
            width = length(string(clms[i]))
            if width > max_width
                width = max_width
            end
            if width > clm_widths[i]
                clm_widths[i] = width
            end
        end
        if size(rows) != (0, 0)
            for i = 1:size(rows)[1], j = 1:size(rows)[2]
                width = length(string(rows[i, j]))
                if width > max_width
                    width = max_width
                end
                if width > clm_widths[j]
                    clm_widths[j] = width
                end
            end
        else
            rows = Matrix(undef, 0, length(clms))
        end
        clm_widths = clm_widths .+ 2
        new(name, clms, rows, clm_widths, max_width)
    end
end

function show(io::IO, tbl::Table)
    if (tbl.name != "")
        println("/ $(tbl.name) \\")
    end
    tbl_width = sum(tbl.clm_widths) + (length(tbl.clms) - 1)
    print("+", rpad("", tbl_width, "-"), "+", "\n|")
    for (i, clm) in enumerate(tbl.clms)
        print(rpad(" $clm ", tbl.clm_widths[i]), "|")
    end
    print("\n+", rpad("", tbl_width, "="), "+");
    for i = 1:size(tbl.rows)[1]
        i == 1 && print("\n|")
        for j = 1:size(tbl.rows)[2]
            print(rpad(" $(tbl.rows[i, j]) ", tbl.clm_widths[j]), "|")
        end
        i != size(tbl.rows)[1] && print("\n|")
    end
    print("\n+", rpad("", tbl_width, "-"), "+")
end

function append!(tbl::Table, row::Vector)
    tbl.rows = [tbl.rows; reshape(row, (1, length(row)))]
    for i = 1:length(row)
        width = length(string(row[i])) + 2
        if width > tbl.max_width
            width = tbl.max_width
        end
        if width > tbl.clm_widths[i]
            tbl.clm_widths[i] = width
        end
    end
    return nothing
end

function insert!(tbl::Table, index::Integer, row::Vector)
    tbl.rows = [tbl.rows[1:(index - 1), :]; reshape(row, (1, length(row))); tbl.rows[index:end, :]]
    for i = 1:length(row)
        width = length(string(row[i])) + 2
        if width > tbl.max_width
            width = tbl.max_width
        end
        if width > tbl.clm_widths[i]
            tbl.clm_widths[i] = width
        end
    end
    return nothing
end

function remove!(tbl::Table, index::Integer)
    row = tbl.rows[index, :]
    tbl.rows = [tbl.rows[1:(index - 1), :]; tbl.rows[(index + 1):end, :]]
    return row
end

end
