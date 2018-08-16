import React from 'react'
import PropTypes from 'prop-types'

import DispatchRequestButton from './DispatchRequestButton'
import { IconButtonContent } from './UpdateButton'

const DeleteButton = ({ initialValues, onSubmit, buttonText, size, ...props }) =>
  <DispatchRequestButton
    buttonContent={<IconButtonContent editIconName="trash" buttonText={buttonText} size={size} />}
    onSubmit={() => onSubmit({ ...initialValues, delete: true })}
    {...props}
  />

DeleteButton.propTypes = {
  onSubmit: PropTypes.func,
  confirmDialog: PropTypes.oneOfType([PropTypes.string, PropTypes.node]),
  initialValues: PropTypes.object,
  buttonText: PropTypes.string,
  size: PropTypes.string,
}

export default DeleteButton
