from pressure import press

def f1(mdot, v2, p2, area, press1, alt, propellantmass):
    press1 = press(alt)
    F = mdot * v2 + (p2 - press1) * area
    if propellantmass < .1:
        thrust = 0
        mDot = 0
    else:
        thrust = F
        #Redundant,but explicit
        mdot = mdot
    return thrust, mDot
