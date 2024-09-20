import bpy
import numpy as np
from math import pi
from mathutils import Vector

# --- Step 1: Load Interface Coordinates ---
path = '/destination_path/'
interface_file = f'{path}/interface-0.95.dat'

# Read the interface data
interface_data = np.loadtxt(interface_file)

# Reshape data: 2 rows correspond to 1 line segment
n_segments = interface_data.shape[0] // 2
segments = interface_data.reshape(n_segments, 2, 2)

# Prepare arrays for x and y coordinates
x_coords = []
y_coords = []

# Loop through the segments and insert NaN between non-continuous segments
for i in range(n_segments):
    # Get the first point of the current segment
    x1, y1 = segments[i, 0, :]
    x2, y2 = segments[i, 1, :]
   
    # Append current segment points to the list
    x_coords.append([y1, y2])  # Use y1 and y2 as x_coords since x is 2nd column
    y_coords.append([x1, x2])  # Use x1 and x2 as y_coords since y is 1st column
   
    # If this is not the last segment, check for discontinuity
    if i < n_segments - 1:
        next_x1, next_y1 = segments[i + 1, 0, :]
       
        # If the next segment does not start where the current one ends, add NaNs
        if not (x2 == next_x1 and y2 == next_y1):
            x_coords.append([np.nan, np.nan])
            y_coords.append([np.nan, np.nan])

# Convert lists into numpy arrays
x_coords = np.concatenate(x_coords)
y_coords = np.concatenate(y_coords)

# --- Step 2: Revolve Interface to Create 3D Mesh ---
def revolve_points_around_axis(x, y, steps=64):
    """ Function to revolve points around the Z-axis to create a 3D surface """
    verts = []
    faces = []
    n_points = len(x)

    # Revolve from 0 to π radians (180 degrees)
    for step in range(steps):
        angle = (pi * step) / (steps - 1)  # Only revolve for 180 degrees (π radians)
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)

        for i in range(n_points):
            verts.append(Vector((x[i] * cos_angle, x[i] * sin_angle, y[i])))

    for step in range(steps - 1):
        for i in range(n_points - 1):
            current = step * n_points + i
            next = (step + 1) * n_points + i
            faces.append([current, current + 1, next + 1, next])

    return verts, faces

# Get 3D vertices and faces
verts, faces = revolve_points_around_axis(x_coords, y_coords)

# Create a new mesh and object in Blender
mesh = bpy.data.meshes.new(name='Revolved_Interface')
mesh.from_pydata(verts, [], faces)
mesh.update()

obj = bpy.data.objects.new('Revolved_Interface', mesh)
bpy.context.collection.objects.link(obj)

# --- Step 3: Create Material with Glossy Effect and Color Mapping ---

# Create a new material
mat = bpy.data.materials.new(name="Glossy_Surface_Material")
mat.use_nodes = True

# Get the nodes for the material
nodes = mat.node_tree.nodes
bsdf = nodes.get("Principled BSDF")

# Add a 'Separate XYZ' node to get the Y-coordinate of each vertex
separate_xyz = nodes.new(type="ShaderNodeSeparateXYZ")

# Add a 'Color Ramp' node to map the Y-coordinate to a color gradient (blue to red)
color_ramp = nodes.new(type="ShaderNodeValToRGB")
color_ramp.color_ramp.interpolation = 'LINEAR'

# Set the colors in the color ramp (blue for -1, red for 1)
color_ramp.color_ramp.elements[0].position = 0.0
color_ramp.color_ramp.elements[0].color = (8./255., 65./255., 123./255., 1)  # Blue for y = -1
color_ramp.color_ramp.elements[1].position = 1.0
color_ramp.color_ramp.elements[1].color = (1, 0, 0, 1)  # Red for y = 1

# Add a 'Mapping' node to normalize the Y-coordinates to range from -1 to 1
mapping_node = nodes.new(type="ShaderNodeMapping")
mapping_node.inputs["Location"].default_value[1] = 0  # No offset on Y-axis
mapping_node.inputs["Scale"].default_value[1] = 2.0 / (max(y_coords) - min(y_coords))  # Normalize to [-1, 1]

# Link the nodes together
mat.node_tree.links.new(separate_xyz.outputs['Y'], mapping_node.inputs['Vector'])
mat.node_tree.links.new(mapping_node.outputs['Vector'], color_ramp.inputs['Fac'])
mat.node_tree.links.new(color_ramp.outputs['Color'], bsdf.inputs['Base Color'])

# Make the surface glossy by setting roughness close to 0
bsdf.inputs['Roughness'].default_value = 0.05  # Glossy effect

# Assign the material to the revolved object
obj.data.materials.append(mat)

# --- Step 5: Set Up Rendering ---

# Set up cycles rendering and enable ray tracing
bpy.context.scene.render.engine = 'CYCLES'
bpy.context.scene.cycles.samples = 128  # Set the number of samples for ray tracing
bpy.context.scene.render.image_settings.file_format = 'PNG'

# Add lighting and camera for ray tracing effects
light_data = bpy.data.lights.new(name="Light", type='POINT')
light_object = bpy.data.objects.new("Light", light_data)
bpy.context.collection.objects.link(light_object)
light_object.location = (0, -5, 5)
light_data.energy = 1000

# Position the camera
camera = bpy.data.objects['Camera']
camera.location = (-4, -6, 1.5)
camera.rotation_euler = (1.39626, 0, -0.523599)

# --- Step 6: Render the Image ---
output_filepath = f"{path}/blender_render_glossy.png"
bpy.context.scene.render.filepath = output_filepath
bpy.ops.render.render(write_still=True)

print(f"Rendered image saved at {output_filepath}")

################################################################################################################################################

import bpy
import numpy as np
from mathutils import Vector
from scipy.interpolate import griddata

# --- Step 1: Load the Data from File ---
path = '/destination_path/'
data_file = f'{path}/blender_data.txt'

# Load the data: x, z, magnitude (optionally r, g, b if precomputed)
data = np.loadtxt(data_file)

# Split into components
x = -1.*data[:, 1]
z = data[:, 0]
magnitudes = data[:, 2]  # If you exported magnitudes
if data.shape[1] > 3:
    r = data[:, 3]
    g = data[:, 4]
    b = data[:, 5]
    a = data[:, 6]

# --- Step 2: Create Mesh from X, Z Coordinates ---
def create_mesh_from_data(x, z):
    verts = []
    faces = []
    grid_size = int(np.sqrt(len(x)))  # Assuming the grid is square

    # Create vertices
    for i in range(len(x)):
        verts.append(Vector((x[i], 0, z[i])))  # Y-axis is 0 for 2D plot

    # Create faces for the grid (quads)
    for i in range(grid_size - 1):
        for j in range(grid_size - 1):
            index = i * grid_size + j
            v1 = index
            v2 = index + 1
            v3 = index + grid_size + 1
            v4 = index + grid_size
            faces.append([v1, v2, v3, v4])

    return verts, faces

verts, faces = create_mesh_from_data(x, z)

# Create the mesh object in Blender
mesh = bpy.data.meshes.new(name='2D_Plane_Mesh')
mesh.from_pydata(verts, [], faces)
mesh.update()

# Create object and link to scene
obj = bpy.data.objects.new('2D_Plane_Object', mesh)
bpy.context.collection.objects.link(obj)

# --- Step 3: Apply Color Mapping ---
# If you have precomputed RGB values in Python, you can use them directly
if data.shape[1] > 3:
    mesh.vertex_colors.new(name="ColorMap")
    color_layer = mesh.vertex_colors["ColorMap"]
   
    for poly in mesh.polygons:
        for loop_index in poly.loop_indices:
            vertex_index = mesh.loops[loop_index].vertex_index
            r_val = r[vertex_index]
            g_val = g[vertex_index]
            b_val = b[vertex_index]
            a_val = a[vertex_index]
            color_layer.data[loop_index].color = (r_val, g_val, b_val, a_val)
else:
    # If not, normalize magnitudes and apply a basic heatmap in Blender
    mesh.vertex_colors.new(name="MagnitudeColors")
    color_layer = mesh.vertex_colors["MagnitudeColors"]

    tp_min = np.min(magnitudes)
    tp_max = np.max(magnitudes)
    norm_magnitudes = (magnitudes - tp_min) / (tp_max - tp_min)

    for poly in mesh.polygons:
        for loop_index in poly.loop_indices:
            vertex_index = mesh.loops[loop_index].vertex_index
            magnitude = norm_magnitudes[vertex_index]
            color_layer.data[loop_index].color = (magnitude, 0, 1 - magnitude, 1.0)  # Example heatmap (blue to red)

# --- Step 4: Set Up Materials ---
mat = bpy.data.materials.new(name="2D_Plane_Material")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]

# Use vertex colors in the material
color_input = mat.node_tree.nodes.new(type="ShaderNodeVertexColor")
color_input.layer_name = "ColorMap" if data.shape[1] > 3 else "MagnitudeColors"

# Link vertex colors to the base color of the material
mat.node_tree.links.new(color_input.outputs["Color"], bsdf.inputs["Base Color"])
mat.node_tree.links.new(color_input.outputs["Alpha"], bsdf.inputs["Alpha"])

mat.blend_method = 'BLEND'
mat.shadow_method = 'NONE'

# Assign the material to the object
obj.data.materials.append(mat)

# --- Step 5: Set Up Rendering ---
# (Same rendering setup as in your original script)
bpy.context.scene.render.engine = 'CYCLES'
bpy.context.scene.cycles.samples = 128  # Set the number of samples for ray tracing
bpy.context.scene.render.image_settings.file_format = 'PNG'

# Position the camera and lighting (adjust as needed)
light_data = bpy.data.lights.new(name="Light", type='POINT')
light_object = bpy.data.objects.new("Light", light_data)
bpy.context.collection.objects.link(light_object)
light_object.location = (0, -5, 5)
light_data.energy = 1000

camera = bpy.data.objects['Camera']
camera.location = (-4, -6, 1.5)
camera.rotation_euler = (1.39626, 0, -0.523599)

# Render the image
output_filepath = f"{path}/blender_render_plane.png"
bpy.context.scene.render.filepath = output_filepath
bpy.ops.render.render(write_still=True)

print(f"Rendered image saved at {output_filepath}")

########################################################################################################################################################

import bpy
import numpy as np
from mathutils import Vector
from scipy.interpolate import griddata

# --- Step 1: Load the Data from File ---
data_file = f'{path}/blender_index_data.txt'

# Load the data: x, z, magnitude (optionally r, g, b if precomputed)
data = np.loadtxt(data_file)

# Split into components
x = data[:, 1]
z = data[:, 0]
magnitudes = data[:, 2]  # If you exported magnitudes
if data.shape[1] > 3:
    r = data[:, 3]
    g = data[:, 4]
    b = data[:, 5]
    a = data[:, 6]

# --- Step 2: Create Mesh from X, Z Coordinates ---
def create_mesh_from_data(x, z):
    verts = []
    faces = []
    grid_size = int(np.sqrt(len(x)))  # Assuming the grid is square

    # Create vertices
    for i in range(len(x)):
        verts.append(Vector((x[i], 0, z[i])))  # Y-axis is 0 for 2D plot

    # Create faces for the grid (quads)
    for i in range(grid_size - 1):
        for j in range(grid_size - 1):
            index = i * grid_size + j
            v1 = index
            v2 = index + 1
            v3 = index + grid_size + 1
            v4 = index + grid_size
            faces.append([v1, v2, v3, v4])

    return verts, faces

verts, faces = create_mesh_from_data(x, z)

# Create the mesh object in Blender
mesh = bpy.data.meshes.new(name='2D_Plane_Mesh_index')
mesh.from_pydata(verts, [], faces)
mesh.update()

# Create object and link to scene
obj = bpy.data.objects.new('2D_Plane_Object_index', mesh)
bpy.context.collection.objects.link(obj)

# --- Step 3: Apply Color Mapping ---
# If you have precomputed RGB values in Python, you can use them directly
if data.shape[1] > 3:
    mesh.vertex_colors.new(name="ColorMap_index")
    color_layer = mesh.vertex_colors["ColorMap_index"]
   
    for poly in mesh.polygons:
        for loop_index in poly.loop_indices:
            vertex_index = mesh.loops[loop_index].vertex_index
            r_val = r[vertex_index]
            g_val = g[vertex_index]
            b_val = b[vertex_index]
            a_val = a[vertex_index]
            color_layer.data[loop_index].color = (r_val, g_val, b_val, a_val)
else:
    # If not, normalize magnitudes and apply a basic heatmap in Blender
    mesh.vertex_colors.new(name="MagnitudeColors")
    color_layer = mesh.vertex_colors["MagnitudeColors"]

    tp_min = np.min(magnitudes)
    tp_max = np.max(magnitudes)
    norm_magnitudes = (magnitudes - tp_min) / (tp_max - tp_min)

    for poly in mesh.polygons:
        for loop_index in poly.loop_indices:
            vertex_index = mesh.loops[loop_index].vertex_index
            magnitude = norm_magnitudes[vertex_index]
            color_layer.data[loop_index].color = (magnitude, 0, 1 - magnitude, 1.0)  # Example heatmap (blue to red)

# --- Step 4: Set Up Materials ---
mat = bpy.data.materials.new(name="2D_Plane_Material_index")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]

# Use vertex colors in the material
color_input = mat.node_tree.nodes.new(type="ShaderNodeVertexColor")
color_input.layer_name = "ColorMap_index" if data.shape[1] > 3 else "MagnitudeColors"

# Link vertex colors to the base color of the material
mat.node_tree.links.new(color_input.outputs["Color"], bsdf.inputs["Base Color"])
mat.node_tree.links.new(color_input.outputs["Alpha"], bsdf.inputs["Alpha"])

mat.blend_method = 'BLEND'
mat.shadow_method = 'NONE'

# Assign the material to the object
obj.data.materials.append(mat)

# --- Step 5: Set Up Rendering ---
# (Same rendering setup as in your original script)
bpy.context.scene.render.engine = 'CYCLES'
bpy.context.scene.cycles.samples = 128  # Set the number of samples for ray tracing
bpy.context.scene.render.image_settings.file_format = 'PNG'

# Position the camera and lighting (adjust as needed)
light_data = bpy.data.lights.new(name="Light", type='POINT')
light_object = bpy.data.objects.new("Light", light_data)
bpy.context.collection.objects.link(light_object)
light_object.location = (0, -5, 5)
light_data.energy = 1000

camera = bpy.data.objects['Camera']
camera.location = (-4, -6, 1.5)
camera.rotation_euler = (1.39626, 0, -0.523599)

# Render the image
output_filepath = f"{path}/blender_render_plane.png"
bpy.context.scene.render.filepath = output_filepath
bpy.ops.render.render(write_still=True)

print(f"Rendered image saved at {output_filepath}")

#########################################################################################
import bpy
import numpy as np
from math import pi
from mathutils import Vector

world = bpy.data.worlds.new(name="world")
bpy.context.scene.world = world
world.color = (0, 0, 0)

# --- Step 5: Set Up Rendering ---
# Set the physical size and DPI
physical_width_cm = 2.4
physical_height_cm = 2.0
dpi = 1200

# Convert physical dimensions from cm to inches (1 cm = 0.393701 inches)
width_inch = physical_width_cm * 0.393701
height_inch = physical_height_cm * 0.393701

# Calculate pixel dimensions based on DPI
pixel_width = int(width_inch * dpi)
pixel_height = int(height_inch * dpi)

# Set the render resolution
bpy.context.scene.render.resolution_x = pixel_width
bpy.context.scene.render.resolution_y = pixel_height
bpy.context.scene.render.resolution_percentage = 100  # Keep at 100%

# Set the render file format to JPEG
bpy.context.scene.render.image_settings.file_format = 'JPEG'
bpy.context.scene.render.image_settings.quality = 100  # Set quality (0-100)

# Position the camera and lighting (adjust as needed)
light_data = bpy.data.lights.new(name="Light", type='POINT')
light_object = bpy.data.objects.new("Light", light_data)
bpy.context.collection.objects.link(light_object)
light_object.location = (0, -5, 5)
light_data.energy = 1000

camera = bpy.data.objects['Camera']
camera.location = (-4, -6, 1.5)
camera.rotation_euler = (1.39626, 0, -0.523599)

# Render the image
output_filepath = f"{path}/blender_render_interpolated_plane_final.jpg"
bpy.context.scene.render.filepath = output_filepath
bpy.ops.render.render(write_still=True)

print(f"Rendered image saved at {output_filepath}")
