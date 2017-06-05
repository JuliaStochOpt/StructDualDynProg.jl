function testniter_3stages(niter, K, maxncuts, cutmode, detectlb)
    if cutmode == :MultiCut
        if K == -1
            if maxncuts == -1
                if detectlb
                    # == 15 on Mac OS, == 17 on Linux
                    @test 15 <= niter <= 17
                else
                    @test niter == 19
                end
            else
                if detectlb
                    # down to 15 and up to 19 on Mac OS, 17 <= ... <= 18 on Linux
                    @test 15 <= niter <= 19
                else
                    @test 19 <= niter <= 21
                end
            end
        else
            if detectlb
                # up to 15 on Mac OS, <= 14 on Linux
                @test 12 <= niter <= 15
            else
                # up to 19 on Mac OS, <= 18 on Linux
                @test 15 <= niter <= 19
            end
        end
    else
        if maxncuts == -1
            if K == -1
                # == 19 on Mac OS, == 21 on Linux
                @test 19 <= niter <= 21
            else
                if detectlb
                    # == 14 on Mac OS, == 16 on Linux
                    @test 14 <= niter <= 16
                else
                    # == 14 on Mac OS, == 17 on Linux
                    @test 14 <= niter <= 17
                end
            end
        else
            if K == -1
                if detectlb
                    # down to 15 Mac OS, >= 19 on Linux
                    @test 15 <= niter <= 24
                else
                    @test 20 <= niter <= 24
                end
            else
                if detectlb
                    # down to 13 Mac OS, >= 15 on Linux
                    @test 13 <= niter <= 19
                else
                    # down to 14 Mac OS, >= 17 on Linux
                    @test 14 <= niter <= 19
                end
            end
        end
    end
end
