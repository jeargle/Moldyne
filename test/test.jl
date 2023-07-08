# John Eargle (mailto: jeargle at gmail.com)
# 2017-2018
# Moldyne test

using moldyne


function print_test_header(test_name)
    border = repeat("*", length(test_name) + 4)
    println(border)
    println("* ", test_name, " *")
    println(border)
end


function test_structure()
    print_test_header("Structure")

    structure1 = Structure("structure1", 2, "test2.pdb")
    structure2 = Structure("structure2", 3, "test2.pdb")

    println(structure1)
    println(structure2)
    println()
end


function main()
    test_structure()
end

main()
