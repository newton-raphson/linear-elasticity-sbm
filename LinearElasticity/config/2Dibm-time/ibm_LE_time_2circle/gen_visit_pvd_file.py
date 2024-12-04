import os
import sys
import glob
startpath = os.path.dirname(os.path.realpath(__file__))
folders = []

# print(sys.argv)
# print(startpath)

res = []
geom = []
for root, dirs, files in os.walk(startpath):
  for f in [f for f in files if f.endswith(".pvtu")]:
    res.append((root, f))
  for f in [f for f in files if f.endswith(".vtp")]:
    geom.append((root, f))

    # print(root, f)
res = sorted(res, key=lambda x: x[1])
geom = sorted(geom, key=lambda x: (os.path.basename(x[0]), x[1]))


ns_res = [r for r in res if "ns" in r[1]]
le_res = [r for r in res if "le" in r[1]]
# if "all" arguments passed
if len(sys.argv) == 2 and sys.argv[1] == "all":
  pass
else:
  # no arguments, only this directory
  ns_res = [r for r in res if "ns" in r[1] and os.path.dirname(r[0]) == startpath]
  le_res = [r for r in res if "le" in r[1] and os.path.dirname(r[0]) == startpath]
  geom = [r for r in geom if os.path.dirname(r[0]) == startpath]

def write_visit_file(name, res):
  if len(res) > 0:
    with open(os.path.join(startpath, name), "w") as f:
      f.write('!NBLOCKS 1\n')
      for r in res:
        f.write("{file}\n".format(time=float(r[1][-10:-4]), file=os.path.join(r[0], r[1])))
        # f.write("!TIME {time}\n{file}\n".format(time=float(r[1][-10:-4]), file=os.path.join(r[0], r[1])))


def write_geom_file(name, geom):
  if len(geom) > 0:
    folders = list(set([g[0] for g in geom]))
    num_geo_per_folder = [len([g[1] for g in geom if g[0] == f]) for f in folders]
    if len(set(num_geo_per_folder)) == 1:
      print("all geom have the same number")
      with open(os.path.join(startpath, name), "w") as f:
        f.write("!NBLOCKS {}\n".format(num_geo_per_folder[0]))
        for i, g in enumerate(geom):
          # if i % num_geo_per_folder[0] == 0:
          #   f.write("!TIME {time}\n".format(time=float(g[1][-9:-3])))
          f.write("{file}\n".format(file=os.path.join(g[0], g[1])))

    else:
      print("some folder have more or less geom files!")
      for no_f, f in zip(num_geo_per_folder, folders):
        print("{} has {} files.".format(f, no_f))


def write_pvd_file(name, res):
  def res_to_ts(res):
    return int(res.replace(".pvtu", "").split("_")[-1])
  if len(res) > 0:
    with open(os.path.join(startpath, name), "w") as f:
      f.write(r"""<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1"
         byte_order="LittleEndian"
         compressor="vtkZLibDataCompressor">
<Collection>
""")
      for r in res:
        f.write("<DataSet timestep=\"{ts}\" group=\"\" part=\"0\" file=\"{f}\"/>\n".format(ts=res_to_ts(r[1]), f=os.path.join(r[0], r[1])))
      f.write(r"""
  </Collection>
</VTKFile>
""")


def write_pvd_file_geo(name, geom):
  def geo_to_ts(geo):
    return int(geo.replace(".vtp", "").split("_")[-1])
  # folders = list(set([g[0] for g in geom]))
  # num_geo_per_folder = [len([g[1] for g in geom if g[0] == f]) for f in folders]
  with open(os.path.join(startpath, name), "w") as f:
    f.write(r"""<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1"
       byte_order="LittleEndian"
       compressor="vtkZLibDataCompressor">
<Collection>
""")
    for i, g in enumerate(geom):
      f.write("<DataSet timestep=\"{ts}\" group=\"\" part=\"0\" file=\"{f}\"/>\n".format(ts=geo_to_ts(g[1]), f=os.path.join(g[0], g[1])))
    # for r in res:
    #   f.write("<DataSet timestep=\"{ts}\" group=\"\" part=\"0\" file=\"{f}\"/>\n".format(ts=res_to_ts(r[1]), f=os.path.join(r[0], r[1])))

    f.write(r"""
  </Collection>
</VTKFile>
""")


write_visit_file("visit_ns.visit", ns_res)
write_visit_file("visit_le.visit", le_res)
write_geom_file("visit_geo.visit", geom)

write_pvd_file("paraview_ns.pvd", ns_res)
write_pvd_file("paraview_le.pvd", le_res)
write_pvd_file_geo("paraview_geo.pvd", geom)
