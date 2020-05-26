import enum

class Components_names(enum.Enum):
    Bulge = 'Bulge'
    Thin_disk = 'Thin disk'
    Thick_disk = 'Thick disk'
    Exponential_disk = 'Exponential disk'
    Dark_halo = 'Dark halo'
    Burkert_halo = 'Burkert halo'
    

def input_component(component):
    
    component_mass = 10
    component_scale_a = 2
    component_scale_b = 2
    print('Set the guess parameters for', component.value)
    try:
        component_mass = float(input('Mass (in M_sun):'))
    except:
        print('No valid Mass for', component, '. It will be taken the default mass:', component_mass, 'M_sun')
    
    
    try:
        component_scale_a = float(input('Radial Scale Length (in kpc):'))
    except:
        print('No valid Radial Scale Length for', component, '. It will be taken the default Radial Scale Lenght:', component_scale_a, 'kpc')

    if component not in [Components_names.Exponential_disk, Components_names.Dark_halo, Components_names.Burkert_halo]:
        try:
            component_scale_b = float(input('Vertical Scale Length (in kpc):'))
        except:
            print('No valid Vertical Scale Length for', component, '. It will be taken the default Vertical Scale Lenght:', component_scale_b, 'kpc')
    
    return component_mass, component_scale_a, component_scale_b
    
    
    
m, a, b = input_component(Components_names.Bulge)
print(m, a, b)

m, a, b = input_component(Components_names.Exponential_disk)
print(m, a, b)