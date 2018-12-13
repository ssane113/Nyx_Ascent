rm -rf CMakeFiles/ CMakeCache.txt

cmake -DASCENT_DIR=/home/users/ssane/Alpine/Ascent/install \
      -DCONDUIT_DIR=/home/users/ssane/Alpine/Conduit/install \
      -DVTKM_DIR=/home/users/ssane/Alpine/VTKM/install \
      -DVTKH_DIR=/home/users/ssane/Alpine/VTKH/install \
      -DVTK_DIR=/home/users/ssane/VTK-6.3.0 \
      .
