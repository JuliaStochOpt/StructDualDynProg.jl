function testniter_3stages(niter, K, maxncuts, cutmode, detectlb)
    if cutmode == :MultiCut
        if K == -1
            if maxncuts == -1
                if detectlb
                    @test niter == 17
                else
                    @test niter == 19
                end
            else
                if detectlb
                    @test 17 <= niter <= 18
                else
                    @test 19 <= niter <= 21
                end
            end
        else
            if detectlb
                @test 12 <= niter <= 14
            else
                @test 15 <= niter <= 18
            end
        end
    else
        if maxncuts == -1
            if K == -1
                @test niter == 21
            else
                if detectlb
                    @test niter == 16
                else
                    @test niter == 17
                end
            end
        else
            if K == -1
                if detectlb
                    @test 19 <= niter <= 24
                else
                    @test 20 <= niter <= 24
                end
            else
                if detectlb
                    @test 15 <= niter <= 19
                else
                    @test 17 <= niter <= 19
                end
            end
        end
    end
end
