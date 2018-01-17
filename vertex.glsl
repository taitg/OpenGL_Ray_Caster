
#version 410

// location indices for these attributes correspond to those specified in the
// InitializeGeometry() function of the main program
layout(location = 0) in vec2 VertexPosition;
layout(location = 1) in vec3 VertexColour;

// output to be interpolated between vertices and passed to the fragment stage
out vec3 Colour;

void main()
{
    // assign vertex position without modification
    gl_Position = vec4(VertexPosition, 0.0, 1.0); // , 0.0

	//float mod = (VertexPosition.z + 2.0) / 3.0;

    // assign output colour to be interpolated
    Colour = VertexColour;
}
