
import xml.etree.ElementTree as ET

# Load the GDML file
tree = ET.parse('protodune_geometry_v6.gdml')
root = tree.getroot()

def return_position(module):
    if module >= 16:
        position_char = "D"
        module -= 16
    else:
        position_char = "U"

    module += 1

    x = y = z = None

    for physvol in root.iter('physvol'):
        volumeref = physvol.find('volumeref')
        if volumeref is not None and volumeref.get('ref') == f"volAuxDet_CRTModule_{position_char}{module}":
            position = physvol.find('position')
            if position is not None:
                x = position.get('x')
                y = position.get('y')
                z = position.get('z')
                #print(x, y, z)
            else:
                print("Position tag not found for module")
            break
    else:
        print("volumeref for module not found")

    return x, y, z, position_char

def return_position2(module, strip):

    if module >= 16:
        position_char = "D"
        module -= 16
    else:
        position_char = "U"

    module += 1
    strip += 1

    x_aux = y_aux = z_aux = None

    for physvol in root.iter('physvol'):
        volumeref = physvol.find('volumeref')
        if volumeref is not None and volumeref.get('ref') == f"volAuxDetSensitive_CRTPaddle_{position_char}{module}_{strip}":
            position = physvol.find('position')
            if position is not None:
                x_aux = position.get('x')
                y_aux = position.get('y')
                z_aux = position.get('z')
                #print(x_aux, y_aux, z_aux)
            else:
                print("Position tag not found for strip")
            break
    else:
        print("volumeref for strip not found")

    return  x_aux, y_aux, z_aux
