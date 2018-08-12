#include "interface.h"

namespace EITS {
	void select(Vector3d *_coord, int _mode, int _type, StlMesh *_stl)
	{
		if (!_stl->getIsLoaded())
			return;

		if (_mode != GLUT_ACTIVE_CTRL) {
			_stl->clearSelection( _type );
		}
		if (true) {
			if (_coord[0] != _coord[1]) {
				_stl->selectWithin(_coord, _type);
			}
			else {
				_stl->selectAt(_coord[8], _type);
			}
		}
	}

	void select(Vector3d *_coord, int _mode, int _type, VolumeMesh *_fem)
	{
		if (!_fem->getIsLoaded())
			return;

		if (_mode != GLUT_ACTIVE_CTRL) {
			_fem->clearSelection(_type);
		}
		if (true) {
			if (_coord[0] != _coord[1]) {
				_fem->selectWithin(_coord, _type);
			}
			else {
				_fem->selectAt(_coord[8], _type);
			}
		}
	}
};
