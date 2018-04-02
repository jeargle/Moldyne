# John Eargle (mailto: jeargle at gmail.com)
# 2017
# Moldyne test

using moldyne


function test_structure()
    println("***")
    println("*** Structure")
    println("***")

    structure1 = Structure("structure1", 2)
    structure2 = Structure("structure2", 3)

    println(structure1)
    println(structure2)
    println()

    setPositions(structure1, "test2.pdb")
    setPositions(structure2, "test2.pdb")

    println(structure1)
    println(structure2)
    println()
end


function main()
    test_structure()
end

main()
