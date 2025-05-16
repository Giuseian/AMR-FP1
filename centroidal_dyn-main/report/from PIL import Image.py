from PIL import Image

# Open the images
img1 = Image.open("C:/Users/giuse/OneDrive/Desktop/GITHUB PROJECTS/AMR-FP1/centroidal_dyn-main/plots/Feet_z walking.png")
img2 = Image.open("C:/Users/giuse/OneDrive/Desktop/GITHUB PROJECTS/AMR-FP1/centroidal_dyn-main/plots/CoM Velocity walking.png")
img3 = Image.open("C:/Users/giuse/OneDrive/Desktop/GITHUB PROJECTS/AMR-FP1/centroidal_dyn-main/plots/Angular Momentum walking.png")

# Resize all to the same width and height
width, height = 600, 400  # Adjust as needed
img1 = img1.resize((width, height))
img2 = img2.resize((width, height))
img3 = img3.resize((width, height))

# Save the adjusted images
img1.save("C:/Users/giuse/OneDrive/Desktop/GITHUB PROJECTS/AMR-FP1/centroidal_dyn-main/plots/Feet_z walking.png")
img2.save("C:/Users/giuse/OneDrive/Desktop/GITHUB PROJECTS/AMR-FP1/centroidal_dyn-main/plots/CoM Velocity walking.png")
img3.save("C:/Users/giuse/OneDrive/Desktop/GITHUB PROJECTS/AMR-FP1/centroidal_dyn-main/plots/Angular Momentum walking.png")
