import React from 'react'
import { shallow, configure } from 'enzyme'
import Adapter from 'enzyme-adapter-react-16'
import { FamilyTableComponent } from './FamilyTable'
import { getVisibleFamiliesInSortedOrder } from '../../utils/visibleFamiliesSelector'

import { STATE1 } from '../../fixtures'

configure({ adapter: new Adapter() })

test('shallow-render without crashing', () => {
  /*
    visibleFamilies: PropTypes.array.isRequired,
    familyGuidToIndividuals: PropTypes.object.isRequired,
   */

  const props = {
    visibleFamilies: getVisibleFamiliesInSortedOrder(STATE1),
  }

  shallow(<FamilyTableComponent {...props} />)
})