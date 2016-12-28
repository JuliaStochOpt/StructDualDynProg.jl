#n = 5
n = 4
m = 3
# Investment cost
#ic = [16, 5, 32, 2, 0] * 8760
if n == 5
    ic = [16, 5, 32, 2, 0]
else
    ic = [16, 5, 32, 2]
end
# Fuel cost
if n == 5
    C = [25, 80, 6.5, 160, 1000]
else
    C = [25, 80, 6.5, 160]
end
# Duration
#T = [8760, 7000, 1500]
T = [8760, 7000, 1500] / 8760
# Height
D10x = diff([0, 3919, 7329, 10315])
Dref = diff([0, 7086, 9004, 11169])
D21 = D10x
D22 = Dref
D2 = [D21  D22]
D31 = D10x
D32 = Dref
D33 = D10x
D34 = Dref
p2 = [0.9, 0.1]
p3 = [0.81, 0.09, 0.09, 0.01]
nstages = 2
