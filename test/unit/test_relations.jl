using Test
using Mimosa

@testset "Group relations parsing" begin
    # Create a temporary TSV file
    mktemp() do path, io
        write(io, "motif\tgroup\n")
        write(io, "m1\tA\n")
        write(io, "m2\tA\n")
        write(io, "m3\tB\n")
        write(io, "m4\tC\n")
        close(io)

        relations = parse_group_relations(path)

        # m1 and m2 are in group A, so they are NOT eligible for each other
        # m1 eligible: m3 (B), m4 (C)
        targets = eligible_targets(relations, "m1")
        @test "m3" in targets
        @test "m4" in targets
        @test "m2" ∉ targets
        @test "m1" ∉ targets

        # m3 (B) eligible: m1 (A), m2 (A), m4 (C)
        targets = eligible_targets(relations, "m3")
        @test "m1" in targets
        @test "m2" in targets
        @test "m4" in targets

        # Groups mapping
        @test relations.groups["m1"] == "A"
        @test relations.groups["m3"] == "B"
    end

    # CSV delimiter detection
    mktemp() do path, io
        write(io, "motif,group\n")
        write(io, "x,G1\n")
        write(io, "y,G2\n")
        close(io)

        relations = parse_group_relations(path)
        targets = eligible_targets(relations, "x")
        @test "y" in targets
        @test isempty(eligible_targets(relations, "y") |> x -> filter(==("x"), x)) == false
    end

    # Custom column names
    mktemp() do path, io
        write(io, "name\tcluster\n")
        write(io, "a\t1\n")
        write(io, "b\t2\n")
        close(io)

        relations = parse_group_relations(path; name_column="name", group_column="cluster")
        targets = eligible_targets(relations, "a")
        @test "b" in targets
    end

    # Known names validation
    mktemp() do path, io
        write(io, "motif\tgroup\n")
        write(io, "known\tA\n")
        write(io, "unknown\tB\n")
        close(io)

        @test_throws ArgumentError parse_group_relations(
            path, known_names=Set(["known", "other"])
        )

        # ignore_missing = true should not throw
        relations = parse_group_relations(
            path, known_names=Set(["known", "other"]), ignore_missing=true
        )
        @test "known" in keys(relations.eligible)
    end

    # Missing required columns
    mktemp() do path, io
        write(io, "motif\tother\n")
        write(io, "m1\tA\n")
        close(io)

        @test_throws ArgumentError parse_group_relations(path)
    end
end
