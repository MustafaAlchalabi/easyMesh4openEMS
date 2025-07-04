# easyMesh4openEMS

## For use
Clone repository with
```
git clone https://github.com/MustafaAlchalabi/easyMesh4openEMS.git
```
to `<path of cloned repository>`.
In order to import the mesher module you have to tell python where to search for it.
For this purpose (as there is no module installation so far) add
```python
sys.path.append('/<path with cloned repository>/easyMesh4openEMS')
import easyMesher
```
to your python script or
```bash
export PYTHONPATH=/<path with cloned repository>/easyMesh4openEMS:$PYTHONPATH
```
to your `.bashrc` or `.profile`, etc.

## For Contributions
The commit masseges have to follow the rules defined in https://www.conventionalcommits.org/en/v1.0.0/
