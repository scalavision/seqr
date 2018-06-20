import React from 'react'
import { shallow, configure } from 'enzyme'
import Adapter from 'enzyme-adapter-react-16'
import { getProjectIndividuals } from '../../selectors'
import { CaseReviewStatusDropdownComponent } from './CaseReviewStatusDropdown'

import { STATE1 } from '../../fixtures'

configure({ adapter: new Adapter() })

test('shallow-render without crashing', () => {
  /*
    individual: PropTypes.object.isRequired,
    updateIndividualsByGuid: PropTypes.func.isRequired,
   */

  const props = {
    individual: getProjectIndividuals(STATE1)[0],
    updateIndividualsByGuid: () => {},
  }

  shallow(<CaseReviewStatusDropdownComponent {...props} />)
})