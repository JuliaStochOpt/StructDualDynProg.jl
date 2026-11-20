function testniter_3stages(niter, maxncuts, cutgen, detectlb)
    if isa(cutgen, StructProg.MultiCutGenerator)
        if maxncuts == -1
            @test 14 <= niter <= 19
        elseif detectlb
            # GLPK can take 22 iterations on 32-bit Linux.
            @test 14 <= niter <= (Sys.WORD_SIZE == 32 ? 22 : 21)
        else
            @test 15 <= niter <= 22
        end
    elseif maxncuts == -1
        @test 19 <= niter <= 23
    elseif detectlb
        @test 15 <= niter <= 27
    else
        @test 20 <= niter <= 27
    end
end
